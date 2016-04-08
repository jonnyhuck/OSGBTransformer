package OSGBtransformer;

/**
 * Class containing static methods to transform 2D or 3Dcoordinates from 
 * geographic wgs84 <-> projected osgb36 (airy)
 * 
 * Please be aware that the "2D" transformations are merely 3D assuming that all 
 *  ellipsoid heights are 0. This has been generally found to have a negligible 
 *  impact upon results, but if precision is important, then please stick 
 *  to 3D transformations
 * 
 * For further information relating top transforming coordinates between OSGB36
 *  and WGS84, please refer to the excellent book by the Ordnance Survey, upon 
 *  which this code is based.
 *  http://www.ordnancesurvey.co.uk/oswebsite/gps/docs/A_Guide_to_Coordinate_Systems_in_Great_Britain.pdf
 * 
 * Copyright (C) 2008 Jonny Huck
 * 
 * This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @version 0.1
 * @author Jonny Huck
 */

public class OSGBGtransformer {
    
    /*
     * class variables to describe OSGB36/Airy and WGS84 coordinate eystems
     */

    //wgs84 Radius of the Equator (Semi-major axis a)
    private static double wgs84_a = 6378137;
    //wgs84 Distance along ellipsoid axis between equator and pole (Semi-minor axis b)
    private static double wgs84_b = 6356752.3142;
    //osgb36 Radius of the Equator (Semi-major axis a)
    private static double osgb36_a = 6377563.396;
    //osgb36 Distance along ellipsoid axis between equator and pole (Semi-minor axis b)
    private static double osgb36_b = 6356256.91;
    //osgb36 Easting True Origin
    private static int E0 = 400000;
    //osgb36 Northing True Origin
    private static int N0 = -100000;
    //osgb36 Scale Factor at natural origin (to minimise distortion with distance from central meridian)
    private static double K0 = 0.9996012717;
    //osgb36 Latitude of natural origin
    private static int PHI0 = 49;
    //osgb36 Longitude of natural origin
    private static int LAM0 = -2;
    //Translation parallel to X
    private static double toOSGB36_tx = -446.488;
    //Translation parallel to Y
    private static double toOSGB36_ty = 125.157;
    //Translation parallel to Z
    private static double toOSGB36_tz = -542.06;
    //Scale change
    private static double toOSGB36_s = 20.4894;
    //Rotation parallel to X
    private static double toOSGB36_rx = -0.1502;
    //Rotation parallel to Y
    private static double toOSGB36_ry = -0.247;
    //Rotation parallel to Z
    private static double toOSGB36_rz = -0.8421;
    //As above but reversed (inverted signs)
    private static double toWGS84_tx = 446.488;
    private static double toWGS84_ty = -125.157;
    private static double toWGS84_tz = 542.06;
    private static double toWGS84_s = -20.4894;
    private static double toWGS84_rx = 0.1502;
    private static double toWGS84_ry = 0.247;
    private static double toWGS84_rz = 0.8421;
    
    /*
     * Public transformation functions
     */
    
    /**
     * Transform a wgs86 lat, lng & height to an osgb36 well known text string
     * @param lat
     * @param lng
     * @param height
     * @return String WKT representing easting and northing
     */
    public static String wgs84_to_osgb36_wkt(double lat, double lng, double height) {

        int easting = wgs84_to_osgb36_e(lat, lng, height);
        int northing = wgs84_to_osgb36_n(lat, lng, height);
        return "POINT(" + Integer.toString(easting) + " " + Integer.toString(northing) + ")";

    }   //wgs84_to_osgb36_wkt
    
    /**
     * Transform a wgs86 lat & lng to an osgb36 well known text string
     * @param lat
     * @param lng
     * @return String WKT representing easting and northing
     */
    public static String wgs84_to_osgb36_wkt(double lat, double lng) {

        int easting = wgs84_to_osgb36_e(lat, lng, 0);
        int northing = wgs84_to_osgb36_n(lat, lng, 0);
        return "POINT(" + Integer.toString(easting) + " " + Integer.toString(northing) + ")";

    }   //wgs84_to_osgb36_wkt

    /**
     * Transform WGS84 Geodetic to OSGB Easting
     * @param lat
     * @param lng
     * @param height
     * @return int easting
     */
    public static int wgs84_to_osgb36_e(double lat, double lng, double height) {

        //Convert to WGS84 Geodetic to WGS84 Cartesian
        double x = latLng_h_to_x(lat, lng, height, wgs84_a, wgs84_b);
        double y = latLng_h_to_y(lat, lng, height, wgs84_a, wgs84_b);
        double z = lat_h_to_z(lat, height, wgs84_a, wgs84_b);

        //Transform to OSGB36 Cartesian
        double x2 = helmert_x(x, y, z, toOSGB36_tx, toOSGB36_ry, toOSGB36_rz, toOSGB36_s);
        double y2 = helmert_y(x, y, z, toOSGB36_ty, toOSGB36_rx, toOSGB36_rz, toOSGB36_s);
        double z2 = helmert_z(x, y, z, toOSGB36_tz, toOSGB36_rx, toOSGB36_ry, toOSGB36_s);

        //Convert to OSGB36 Geodetic
        double lat2 = xyz_to_lat(x2, y2, z2, osgb36_a, osgb36_b);
        double lng2 = xyz_to_lng(x2, y2);
        double h2 = xyz_to_h(x2, y2, z2, osgb36_a, osgb36_b);

        //Convert to OSGB36 E (Remove sub-metre data)
        return (int) latLng_to_e(lat2, lng2, osgb36_a, osgb36_b);

    }   //wgs84_to_osgb36_e
    
    /**
     * Transform WGS84 Geodetic to OSGB Easting
     * @param lat
     * @param lng
     * @return int easting
     */
    public static int wgs84_to_osgb36_e(double lat, double lng) {

        //Convert to WGS84 Geodetic to WGS84 Cartesian
        double x = latLng_h_to_x(lat, lng, 0, wgs84_a, wgs84_b);
        double y = latLng_h_to_y(lat, lng, 0, wgs84_a, wgs84_b);
        double z = lat_h_to_z(lat, 0, wgs84_a, wgs84_b);

        //Transform to OSGB36 Cartesian
        double x2 = helmert_x(x, y, z, toOSGB36_tx, toOSGB36_ry, toOSGB36_rz, toOSGB36_s);
        double y2 = helmert_y(x, y, z, toOSGB36_ty, toOSGB36_rx, toOSGB36_rz, toOSGB36_s);
        double z2 = helmert_z(x, y, z, toOSGB36_tz, toOSGB36_rx, toOSGB36_ry, toOSGB36_s);

        //Convert to OSGB36 Geodetic
        double lat2 = xyz_to_lat(x2, y2, z2, osgb36_a, osgb36_b);
        double lng2 = xyz_to_lng(x2, y2);
        double h2 = xyz_to_h(x2, y2, z2, osgb36_a, osgb36_b);

        //Convert to OSGB36 E (Remove sub-metre data)
        return (int) latLng_to_e(lat2, lng2, osgb36_a, osgb36_b);

    }   //wgs84_to_osgb36_e

    /**
     * Transform WGS84 Geodetic to OSGB Northing
     * @param lat
     * @param lng
     * @param height
     * @return int easting converted from the supplied wgs84 coordinates
     */
    public static int wgs84_to_osgb36_n(double lat, double lng, double height) {

        //Convert to WGS84 Geodetic to WGS84 Cartesian
        double x = latLng_h_to_x(lat, lng, height, wgs84_a, wgs84_b);
        double y = latLng_h_to_y(lat, lng, height, wgs84_a, wgs84_b);
        double z = lat_h_to_z(lat, height, wgs84_a, wgs84_b);

        //Transform to OSGB36 Cartesian
        double x2 = helmert_x(x, y, z, toOSGB36_tx, toOSGB36_ry, toOSGB36_rz, toOSGB36_s);
        double y2 = helmert_y(x, y, z, toOSGB36_ty, toOSGB36_rx, toOSGB36_rz, toOSGB36_s);
        double z2 = helmert_z(x, y, z, toOSGB36_tz, toOSGB36_rx, toOSGB36_ry, toOSGB36_s);

        //Convert to OSGB36 Geodetic
        double lat2 = xyz_to_lat(x2, y2, z2, osgb36_a, osgb36_b);
        double lng2 = xyz_to_lng(x2, y2);
        double h2 = xyz_to_h(x2, y2, z2, osgb36_a, osgb36_b);

        //Convert to OSGB36 N (Remove sub-metre data)
        return (int) latLng_to_n(lat2, lng2, osgb36_a, osgb36_b);

    }   //wgs84_to_osgb36_n

    /**
     * Transform WGS84 Geodetic to OSGB Northing
     * @param lat
     * @param lng
     * @return int easting converted from the supplied wgs84 coordinates
     */
    public static int wgs84_to_osgb36_n(double lat, double lng) {

        //Convert to WGS84 Geodetic to WGS84 Cartesian
        double x = latLng_h_to_x(lat, lng, 0, wgs84_a, wgs84_b);
        double y = latLng_h_to_y(lat, lng, 0, wgs84_a, wgs84_b);
        double z = lat_h_to_z(lat, 0, wgs84_a, wgs84_b);

        //Transform to OSGB36 Cartesian
        double x2 = helmert_x(x, y, z, toOSGB36_tx, toOSGB36_ry, toOSGB36_rz, toOSGB36_s);
        double y2 = helmert_y(x, y, z, toOSGB36_ty, toOSGB36_rx, toOSGB36_rz, toOSGB36_s);
        double z2 = helmert_z(x, y, z, toOSGB36_tz, toOSGB36_rx, toOSGB36_ry, toOSGB36_s);

        //Convert to OSGB36 Geodetic
        double lat2 = xyz_to_lat(x2, y2, z2, osgb36_a, osgb36_b);
        double lng2 = xyz_to_lng(x2, y2);
        double h2 = xyz_to_h(x2, y2, z2, osgb36_a, osgb36_b);

        //Convert to OSGB36 N (Remove sub-metre data)
        return (int) latLng_to_n(lat2, lng2, osgb36_a, osgb36_b);

    }   //wgs84_to_osgb36_n

    /**
     * Transform osgb36 x y z to a wgs84 well known text string
     * @param easting
     * @param northing
     * @param z
     * @return String WKT representing wgs84 latitude and longitude
     */
    public static String osgb36_to_wgs84_wkt(double easting, double northing, double z) {

        //Return wkt lat long (wgs84) from OS
        double lat = osgb36_to_wgs84_lat(easting, northing, z);
        double lng = osgb36_to_wgs84_lng(easting, northing, z);
        return "POINT(" + Double.toString(lat) + " " + Double.toString(lng) + ")";

    }   //osgb36_to_wgs84_wkt

    /**
     * Transform osgb36 x y z to a wgs84 well known text string
     * @param easting
     * @param northing
     * @return String WKT representing wgs84 latitude and longitude
     */
    public static String osgb36_to_wgs84_wkt(double easting, double northing) {

        //Return wkt lat long (wgs84) from OS
        double lat = osgb36_to_wgs84_lat(easting, northing, 0);
        double lng = osgb36_to_wgs84_lng(easting, northing, 0);
        return "POINT(" + Double.toString(lat) + " " + Double.toString(lng) + ")";

    }   //osgb36_to_wgs84_wkt

    /**
     * Transform osgb36 to a wgs84 latitude value
     * @param easting
     * @param northing
     * @param z
     * @return double wgs84 latitude
     */
    public static double osgb36_to_wgs84_lat(double easting, double northing, double z) {

        //Convert to OSGB36 Geodetic
        double lat = en_to_lat(easting, northing, osgb36_a, osgb36_b);
        double lng = en_to_lng(easting, northing, osgb36_a, osgb36_b); //Convert to OSGB36 Cartesian
        double x = latLng_h_to_x(lat, lng, z, osgb36_a, osgb36_b);
        double y = latLng_h_to_y(lat, lng, z, osgb36_a, osgb36_b);
        double cz = lat_h_to_z(lat, z, osgb36_a, osgb36_b); //Cartesian z

        //Transform to WGS84 Cartesian
        double x2 = helmert_x(x, y, cz, toWGS84_tx, toWGS84_ry, toWGS84_rz, toWGS84_s);
        double y2 = helmert_y(x, y, cz, toWGS84_ty, toWGS84_rx, toWGS84_rz, toWGS84_s);
        double z2 = helmert_z(x, y, cz, toWGS84_tz, toWGS84_rx, toWGS84_ry, toWGS84_s); //Convert to WGS84 Geodetic

        return xyz_to_lat(x2, y2, z2, wgs84_a, wgs84_b);

    }   //osgb36_to_wgs84_lat

    /**
     * Transform osgb36 to a wgs84 latitude value
     * @param easting
     * @param northing
     * @return double wgs84 latitude
     */
    public static double osgb36_to_wgs84_lat(double easting, double northing) {

        //Convert to OSGB36 Geodetic
        double lat = en_to_lat(easting, northing, osgb36_a, osgb36_b);
        double lng = en_to_lng(easting, northing, osgb36_a, osgb36_b); //Convert to OSGB36 Cartesian
        double x = latLng_h_to_x(lat, lng, 0, osgb36_a, osgb36_b);
        double y = latLng_h_to_y(lat, lng, 0, osgb36_a, osgb36_b);
        double cz = lat_h_to_z(lat, 0, osgb36_a, osgb36_b); //Cartesian z

        //Transform to WGS84 Cartesian
        double x2 = helmert_x(x, y, cz, toWGS84_tx, toWGS84_ry, toWGS84_rz, toWGS84_s);
        double y2 = helmert_y(x, y, cz, toWGS84_ty, toWGS84_rx, toWGS84_rz, toWGS84_s);
        double z2 = helmert_z(x, y, cz, toWGS84_tz, toWGS84_rx, toWGS84_ry, toWGS84_s); //Convert to WGS84 Geodetic

        return xyz_to_lat(x2, y2, z2, wgs84_a, wgs84_b);

    }   //osgb36_to_wgs84_lat

    /**
     * Transforms osgb36 x y z to wgs84 longitude
     * @param easting
     * @param northing
     * @param z
     * @return double wgs84 longitude
     */
    public static double osgb36_to_wgs84_lng(double easting, double northing, double z){

        //Convert to OSGB36 Geodetic
        double lat = en_to_lat(easting, northing, osgb36_a, osgb36_b);
        double lng = en_to_lng(easting, northing, osgb36_a, osgb36_b);

        //Convert to OSGB36 Cartesian
        double x = latLng_h_to_x(lat, lng, z, osgb36_a, osgb36_b);
        double y = latLng_h_to_y(lat, lng, z, osgb36_a, osgb36_b);
        double cz = lat_h_to_z(lat, z, osgb36_a, osgb36_b);  //Cartesian z

        //Transform to WGS84 Cartesian
        double x2 = helmert_x(x, y, cz, toWGS84_tx, toWGS84_ry, toWGS84_rz, toWGS84_s);
        double y2 = helmert_y(x, y, cz, toWGS84_ty, toWGS84_rx, toWGS84_rz, toWGS84_s);
        double z2 = helmert_z(x, y, cz, toWGS84_tz, toWGS84_rx, toWGS84_ry, toWGS84_s);

        //Convert to WGS84 Geodetic
        return xyz_to_lng(x2, y2);

    }   //osgb36_to_wgs84_lng

    /**
     * Transforms osgb36 x y z to wgs84 longitude
     * @param easting
     * @param northing
     * @return double  longitude
     */
    public static double osgb36_to_wgs84_lng(double easting, double northing){

        //Convert to OSGB36 Geodetic
        double lat = en_to_lat(easting, northing, osgb36_a, osgb36_b);
        double lng = en_to_lng(easting, northing, osgb36_a, osgb36_b);

        //Convert to OSGB36 Cartesian
        double x = latLng_h_to_x(lat, lng, 0, osgb36_a, osgb36_b);
        double y = latLng_h_to_y(lat, lng, 0, osgb36_a, osgb36_b);
        double cz = lat_h_to_z(lat, 0, osgb36_a, osgb36_b);  //Cartesian z

        //Transform to WGS84 Cartesian
        double x2 = helmert_x(x, y, cz, toWGS84_tx, toWGS84_ry, toWGS84_rz, toWGS84_s);
        double y2 = helmert_y(x, y, cz, toWGS84_ty, toWGS84_rx, toWGS84_rz, toWGS84_s);
        double z2 = helmert_z(x, y, cz, toWGS84_tz, toWGS84_rx, toWGS84_ry, toWGS84_s);

        //Convert to WGS84 Geodetic
        return xyz_to_lng(x2, y2);

    }   //osgb36_to_wgs84_lng

    /*
     * Conversion (not transformation) functions
     */

    /**
     * Convert a cartesian coordinate pair to lat (not transformation!)
     * @param east
     * @param North
     * @param a
     * @param b
     * @return double latitude
     */
    private static double en_to_lat(double east, double North, double a, double b){

        //PRIMARY ELLIPSOID PARAMETERS
        //Convert origins to radians
        double RadPHI0 = PHI0 * (Math.PI / 180);
        double RadLAM0 = LAM0 * (Math.PI / 180);

        //Apply scale factor to axes
        double aK0 = a * K0;
        double bK0 = b * K0;

        //DERIVED ELLIPSOID PARAMETERS
        //Second eccentricity of source ellipsoid
        double e2 = ((Math.pow(aK0, 2)) - (Math.pow(bK0, 2))) / (Math.pow(aK0, 2));
        //Second flattening
        double n = (aK0 - bK0) / (aK0 + bK0);
        //
        double Et = east - E0;
        //Calculate initial value for latitude
        double PHId = initialLat(North, N0, aK0, RadPHI0, n, bK0);
        //Transverse radius of prime vertical
        double nu = aK0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHId)), 2))));
        //Transverse radius of curvature of meridian
        double rho = (nu * (1 - e2)) / (1 - (e2 * Math.pow(Math.sin(PHId), 2)));
        //
        double eta2 = (nu / rho) - 1;

        //Calculate latitude
        double VII = (Math.tan(PHId)) / (2 * rho * nu);
        double VIII = ((Math.tan(PHId)) / (24 * rho * Math.pow(nu, 3)) * 
                (5 + (3 * Math.pow((Math.tan(PHId)), 2)) + eta2 - (9 * eta2 * 
                Math.pow((Math.tan(PHId)), 2))));
        double IX = ((Math.tan(PHId)) / (720 * rho * Math.pow(nu, 5))) * 
                (61 + (90 * Math.pow((Math.tan(PHId)), 2)) + 
                (45 * Math.pow((Math.tan(PHId)), 4)));

        return ((180 / Math.PI) * (PHId - (Math.pow(Et, 2) * VII) + 
                (Math.pow(Et, 4) * VIII) - (Math.pow(Et, 6) * IX)));

    }   //en_to_lat


    /**
     * Convert a cartesian coordinate pair to lng (not transformation!)
     * @param East
     * @param North
     * @param a
     * @param b
     * @return double longitude
     */
    private static double en_to_lng(double East, double North, double a, double b){

        //PRIMARY ELLIPSOID PARAMETERS
        //Convert origins to radians
        double RadPHI0 = PHI0 * (Math.PI / 180);
        double RadLAM0 = LAM0 * (Math.PI / 180);
        //Apply scale factor to axes
        double aK0 = a * K0;
        double bK0 = b * K0;

        //DERIVED ELLIPSIOD PARAMETERS
        //Second eccentricity of source ellipsoid
        double e2 = (Math.pow(aK0, 2) - Math.pow(bK0, 2)) / Math.pow(aK0, 2);
        //Second flattening
        double n = (aK0 - bK0) / (aK0 + bK0);
        //
        double Et = East - E0;
        //Calculate initial value for latitude
        double PHId = initialLat(North, N0, aK0, RadPHI0, n, bK0);
        //Transverse radius of prime vertical
        double nu = aK0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHId)), 2))));
        //Transverse radius of curvature of meridian
        double rho = (nu * (1 - e2)) / (1 - (e2 * Math.pow((Math.sin(PHId)), 2)));
        //
        double eta2 = (nu / rho) - 1;

        //Calculate longitude
        double X = Math.pow((Math.cos(PHId)), -1) / nu;
        double XI = (Math.pow((Math.cos(PHId)), -1) / (6 * Math.pow(nu, 3))) * 
                ((nu / rho) + (2 * Math.pow((Math.tan(PHId)), 2)));
        double XII = (Math.pow((Math.cos(PHId)), -1) / (120 * Math.pow(nu, 5))) * 
                (5 + (28 * Math.pow((Math.tan(PHId)), 2)) + (24 * Math.pow((Math.tan(PHId)), 4)));
        double XIIA = (Math.pow((Math.cos(PHId)), -1) / (5040 * Math.pow(nu, 7))) * 
                (61 + (662 * Math.pow((Math.tan(PHId)), 2)) + (1320 * Math.pow((Math.tan(PHId)), 4)) + 
                (720 * Math.pow((Math.tan(PHId)), 6)));

        return ((180 / Math.PI) * (RadLAM0 + (Et * X) - (Math.pow(Et, 3) * XI) + 
                (Math.pow(Et, 5) * XII) - (Math.pow(Et, 7) * XIIA)));

    }   //en_to_lng


    /**
     * Convert a lat, lng to a cartesian easting (not transformation!)
     * @param lat
     * @param lng
     * @param a
     * @param b
     * @return double easting
     */
    private static double latLng_to_e(double lat, double lng, double a, double b){

        //PRIMARY ELLIPSOID PARAMETERS
        //Convert lat, lng and origins to rads
        double RadPHI = lat * (Math.PI / 180);
        double RadLAM = lng * (Math.PI / 180);
        double RadPHI0 = PHI0 * (Math.PI / 180);
        double RadLAM0 = LAM0 * (Math.PI / 180);

        //Apply scale factor to axes
        double aK0 = a * K0;
        double bK0 = b * K0;

        //DERIVED ELLIPSIOD PARAMETERS
        //Second eccentricity of source ellipsoid
        double e2 = (Math.pow(aK0, 2) - Math.pow(bK0, 2)) / Math.pow(aK0, 2);
        //Second flattening
        double n = (aK0 - bK0) / (aK0 + bK0);
        //Transverse radius of prime vertical
        double nu = aK0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));
        //Transverse radius of curvature of meridian
        double rho = (nu * (1 - e2)) / (1 - (e2 * Math.pow(Math.sin(RadPHI), 2)));
        //
        double eta2 = (nu / rho) - 1;

        //Calculate easting
        double p = RadLAM - RadLAM0;
        double IV = nu * (Math.cos(RadPHI));
        double V = (nu / 6) * Math.pow((Math.cos(RadPHI)), 3) * ((nu / rho) - 
                (Math.pow(Math.tan(RadPHI), 2)));
        double VI = (nu / 120) * Math.pow((Math.cos(RadPHI)), 5) * 
                (5 - (18 * Math.pow((Math.tan(RadPHI)), 2)) + Math.pow((Math.tan(RadPHI)), 4) + 
                (14 * eta2) - (58 * Math.pow((Math.tan(RadPHI)), 2) * eta2));

        return (E0 + (p * IV) + (Math.pow(p, 3) * V) + (Math.pow(p, 5) * VI));
        
    }   //latLng_to_e


    /**
     * Convert a lat, lng to a cartesian easting
     * @param lat
     * @param lng
     * @param a
     * @param b
     * @return double northing
     */
    private static double latLng_to_n(double lat, double lng, double a, double b){

        //PRIMARY ELLIPSOID PARAMETERS
        //Convert lat, lng and origins to rads
        double RadPHI = lat * (Math.PI / 180);
        double RadLAM = lng * (Math.PI / 180);
        double RadPHI0 = PHI0 * (Math.PI / 180);
        double RadLAM0 = LAM0 * (Math.PI / 180);

        //Apply scale factor to axes
        double aK0 = a * K0;
        double bK0 = b * K0;

        //DERIVED ELLIPSIOD PARAMETERS
        //Second eccentricity of source ellipsoid
        double e2 = (Math.pow(aK0, 2) - Math.pow(bK0, 2)) / Math.pow(aK0, 2);
        //Second flattening
        double n = (aK0 - bK0) / (aK0 + bK0);
        //Transverse radius of prime vertical
        double nu = aK0 / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));
        //Transverse radius of curvature of meridian
        double rho = (nu * (1 - e2)) / (1 - (e2 * Math.pow(Math.sin(RadPHI), 2)));
        //
        double eta2 = (nu / rho) - 1;

        //Calculate northing
        double p = RadLAM - RadLAM0;
        double M = Marc(bK0, n, RadPHI0, RadPHI);
        double I = M + N0;
        double II = (nu / 2) * (Math.sin(RadPHI)) * (Math.cos(RadPHI));
        double III = ((nu / 24) * (Math.sin(RadPHI)) * Math.pow((Math.cos(RadPHI)), 3)) * 
                (5 - Math.pow((Math.tan(RadPHI)), 2) + (9 * eta2));
        double IIIA = ((nu / 720) * (Math.sin(RadPHI)) * Math.pow((Math.cos(RadPHI)), 5)) * 
                (61 - (58 * Math.pow((Math.tan(RadPHI)), 2)) + Math.pow((Math.tan(RadPHI)), 4));

        return (I + (Math.pow(p, 2) * II) + (Math.pow(p, 4) * III) + (Math.pow(p, 6) * IIIA));

    }   //latLng_to_n


    /**
     * Convert height between WGS84 and Airy
     * @param PHI
     * @param LAM
     * @param H
     * @param a
     * @param b
     * @return
     */
    private static double latLng_h_to_x(double PHI, double LAM, double H, double a, double b){

        //PRIMARY ELLIPSOID PARAMETERS
        //Convert angle measures to radians
        double RadPHI = PHI * (Math.PI / 180);
        double RadLAM = LAM * (Math.PI / 180);

        //DERIVED ELLIPSOID PARAMETERS
        //Eccentricity squared
        double e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
        //Transverse radius of prime vertical
        double V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(RadPHI)), 2))));

        //Calculate X
        return (V + H) * (Math.cos(RadPHI)) * (Math.cos(RadLAM));

    }   //latLng_h_to_x


    /**
     * Get northing from 3D lat lng 
     * @param PHI
     * @param LAM
     * @param H
     * @param a
     * @param b
     * @return
     */
    private static double latLng_h_to_y(double PHI, double LAM, double H, double a, double b){

        //PRIMARY ELLIPSOID PARAMETERS
        //Convert angle measures to radians
        double RadPHI = PHI * (Math.PI / 180);
        double RadLAM = LAM * (Math.PI / 180);

        //DERIVED ELLIPSOID PARAMETERS
        //Eccentricity squared
        double e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
        //Transverse radius of prime vertical
        double V = a / (Math.sqrt(1 - (e2 * Math.pow(Math.sin(RadPHI), 2))));

        //Calculate Y
        return (V + H) * (Math.cos(RadPHI)) * (Math.sin(RadLAM));
    }   //latLng_h_to_y


    /**
     * Get lat from cartesian coords
     * @param X
     * @param Y
     * @param Z
     * @param a
     * @param b
     * @return
     */
    private static double xyz_to_lat(double X, double Y, double Z, double a, double b){

        //
        double RootXYSqr = Math.sqrt(Math.pow(X, 2) + Math.pow(Y, 2));

        //DERIVED ELLIPSOID PARAMETERS
        //Eccentricity squared
        double e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
        //
        double PHI1 = Math.atan(Z / (RootXYSqr * (1 - e2)));

        //Calculate Lat (an iterative process)
        double PHI = Iterate_XYZ_to_Lat(a, e2, PHI1, Z, RootXYSqr);

        //Return to degrees and return
        return PHI * (180 / Math.PI);

    }   //xyz_to_lat


    /**
     * Calculates latitude iteratively until the difference between two iterations is
     *  sufficiently precise
     * @param a
     * @param e2
     * @param PHI1
     * @param Z
     * @param RootXYSqr
     * @return double lat
     */
    private static double Iterate_XYZ_to_Lat(double a, double e2, double PHI1, double Z, double RootXYSqr){

        //Transverse radius of prime vertcal
        double V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHI1)), 2))));
        //Calculate initial Lat Value
        double PHI2 = Math.atan((Z + (e2 * V * (Math.sin(PHI1)))) / RootXYSqr);

        //Iterate towards required accuracy
        while (Math.abs(PHI1 - PHI2) > 0.000000001){
            PHI1 = PHI2;
            V = a / (Math.sqrt(1 - (e2 * Math.pow((Math.sin(PHI1)), 2))));
            PHI2 = Math.atan((Z + (e2 * V * (Math.sin(PHI1)))) / RootXYSqr);
        }

        //Return the resulting value
        return PHI2;

    }   //Iterate_XYZ_to_Lat

    
    /**
     * Get lng (quadrant based)
     * @param X
     * @param Y
     * @return double lng
     */
    private static double xyz_to_lng(double X, double Y){

        //Output dependent on equatorial quadrant as determined by the signs of X and Y
        if (X >= 0){  //Longitude is in the W90 thru 0 to E90 hemisphere
            return (Math.atan(Y / X)) * (180 / Math.PI);
        }

        if (X < 0 && Y >= 0){   //Longitude is in the E90 to E180 quadrant
            return ((Math.atan(Y / X)) * (180 / Math.PI)) + 180;
        }

        if (X < 0 && Y < 0){     //Longitude is in the E180 to W90 quadrant
            return ((Math.atan(Y / X)) * (180 / Math.PI)) - 180;
        }

        //On error
        return -1;

    }   //xyz_to_lng


    /*
     * Ellipsoidal Height Transformation Functions
     */

    /**
     * convert ellipsiodal height 
     * @param PHI
     * @param H
     * @param a
     * @param b
     * @return double z
     */
    private static double lat_h_to_z(double PHI, double H, double a, double b){

        //Convert angle measures to radians
        double RadPHI = PHI * (Math.PI / 180);

        //Eccentricity squared
        double e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
        //Transverse radius of prime vertical
        double V = a / (Math.sqrt(1 - (e2 * (Math.pow(Math.sin(RadPHI), 2)))));

        //Z value
        return ((V * (1 - e2)) + H) * (Math.sin(RadPHI));

    }   //lat_h_to_z


    /**
     * Get ellipsiodal height value from 3D cartesian coordinates
     * @param X
     * @param Y
     * @param Z
     * @param a
     * @param b
     * @return double height
     */
    private static double xyz_to_h(double X, double Y, double Z, double a, double b){

        //Calculate Latitude
        double PHI = xyz_to_lat(X, Y, Z, a, b);
        //Convert lat to radians
        double RadPHI = PHI * (Math.PI / 180);

        //
        double RootXYSqr = Math.sqrt(Math.pow(X, 2) + Math.pow(Y, 2));
        //Eccentricity squared
        double e2 = (Math.pow(a, 2) - Math.pow(b, 2)) / Math.pow(a, 2);
        //Transverse radius of prime vertical
        double V = a / (Math.sqrt(1 - (e2 * Math.pow(Math.sin(RadPHI), 2))));

        //Ellipsoidal Height
        return (RootXYSqr / Math.cos(RadPHI)) - V;

    }   //xyz_to_h


    /*
     * Transformation functions
     */
    
    /**
     * Transform a coordinate between two datums
     * @param X
     * @param Y
     * @param Z
     * @param DX
     * @param Y_Rot
     * @param Z_Rot
     * @param s
     * @return double x
     */
    private static double helmert_x(double X, double Y, double Z, double DX, double Y_Rot, double Z_Rot, double s){
        //Calculate Helmert transformed X coordinate.
        //Input: cartesian XYZ coords (X,Y,Z), X translation (DX) all in meters ;
        //       Y and Z rotations in seconds of arc (Y_Rot, Z_Rot) and scale in ppm (s).

        //Convert ppm scale to a factor
        double sFactor = s * 0.000001;
        //Convert rotations to radians
        double RadY_Rot = (Y_Rot / 3600) * (Math.PI / 180);
        double RadZ_Rot = (Z_Rot / 3600) * (Math.PI / 180);

        //Compute transformed X coord
        return X + (X * sFactor) - (Y * RadZ_Rot) + (Z * RadY_Rot) + DX;
    
    }   //helmert_x


    /**
     * Transform a coordinate between two datums
     * @param X
     * @param Y
     * @param Z
     * @param DY
     * @param X_Rot
     * @param Z_Rot
     * @param s
     * @return double y
     */
    private static double helmert_y(double X, double Y, double Z, double DY, double X_Rot, double Z_Rot, double s){
        //Computed Helmert transformed Y coordinate.
        //Input: cartesian XYZ coords (X,Y,Z), Y translation (DY) all in meters ; _
        //       X and Z rotations in seconds of arc (X_Rot, Z_Rot) and scale in ppm (s).

        //Convert ppm scale to a factor
        double sFactor = s * 0.000001;
        //Convert rotations to radians
        double RadX_Rot = (X_Rot / 3600) * (Math.PI / 180);
        double RadZ_Rot = (Z_Rot / 3600) * (Math.PI / 180);

        //Compute transformed Y coord
        return (X * RadZ_Rot) + Y + (Y * sFactor) - (Z * RadX_Rot) + DY;

    }   //helmert_y

    
    /**
     * Transform a coordinate between two datums
     * @param X
     * @param Y
     * @param Z
     * @param DZ
     * @param X_Rot
     * @param Y_Rot
     * @param s
     * @return double z
     */
    private static double helmert_z(double X, double Y, double Z, double DZ, double X_Rot, double Y_Rot, double s){
        //Computed Helmert transformed Z coordinate.

        //Convert ppm scale to a factor
        double sfactor = s * 0.000001;
        //Convert rotations to radians
        double RadX_Rot = (X_Rot / 3600) * (Math.PI / 180);
        double RadY_Rot = (Y_Rot / 3600) * (Math.PI / 180);

        //Compute transformed Z coord
        return (-1 * X * RadY_Rot) + (Y * RadX_Rot) + Z + (Z * sfactor) + DZ;

    }   // helmert_z


    /*
     * Utility functions
     */
    
    /**
     * Iterative function to calculate Latitude
     * @param North
     * @param n0
     * @param afo
     * @param PHI0
     * @param n
     * @param bfo
     * @return double lat
     */
    private static double initialLat(double North, double n0, double afo, double PHI0, double n, double bfo){
        
        //Calculate initial value for Latitude
        double PHI1 = ((North - n0) / afo) + PHI0;
        double M = Marc(bfo, n, PHI0, PHI1);
        double PHI2 = ((North - n0 - M) / afo) + PHI1;

        while (Math.abs(North - n0 - M) > 0.00001){
            PHI2 = ((North - n0 - M) / afo) + PHI1;
            M = Marc(bfo, n, PHI0, PHI2);
            PHI1 = PHI2;
        }
        
        return PHI2;
    }

    
    /**
     * Calculate meridional arc distance from equator to latitude of first standard parallel
     * @param bK0
     * @param n
     * @param PHI0
     * @param PHI
     * @return double
     */
    private static double Marc(double bK0, double n, double PHI0, double PHI){
        //(Elliptic integral) - Spherical formula
        return (bK0 * (((1 + n + ((5 / 4) * Math.pow(n, 2)) + ((5 / 4) * Math.pow(n, 3))) * (PHI - PHI0))
        - (((3 * n) + (3 * Math.pow(n, 2)) + ((21 / 8) * Math.pow(n, 3))) * (Math.sin(PHI - PHI0)) * (Math.cos(PHI + PHI0)))
        + ((((15 / 8) * Math.pow(n, 2)) + ((15 / 8) * Math.pow(n, 3))) * (Math.sin(2 * (PHI - PHI0))) * (Math.cos(2 * (PHI + PHI0))))
        - (((35 / 24) * Math.pow(n, 3)) * (Math.sin(3 * (PHI - PHI0))) * (Math.cos(3 * (PHI + PHI0))))));

    }   //Marc
    
}   //class