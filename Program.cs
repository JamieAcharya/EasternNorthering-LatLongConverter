using System;

namespace CoordinateConverter
{
    /// <summary>
    /// Converts Easting/Northing coordinates to Latitude/Longitude
    /// Supports UTM and British National Grid (BNG) coordinate systems
    /// </summary>
    public class CoordinateConverter
    {
        // WGS84 ellipsoid parameters
        private const double WGS84_A = 6378137.0;              // Semi-major axis (m)
        private const double WGS84_F = 1.0 / 298.257223563;    // Flattening
        private const double WGS84_E2 = 2 * WGS84_F - WGS84_F * WGS84_F; // Eccentricity squared

        // OSGB36 ellipsoid parameters (for British National Grid)
        private const double OSGB36_A = 6377563.396;
        private const double OSGB36_B = 6356256.909;
        private const double OSGB36_F0 = 0.9996012717;
        private const double OSGB36_LAT0 = 49.0 * Math.PI / 180.0;
        private const double OSGB36_LON0 = -2.0 * Math.PI / 180.0;
        private const double OSGB36_N0 = -100000.0;
        private const double OSGB36_E0 = 400000.0;

        public enum CoordinateSystem
        {
            UTM,
            BritishNationalGrid
        }

        public class LatLon
        {
            public double Latitude { get; set; }
            public double Longitude { get; set; }

            public override string ToString()
            {
                return $"Lat: {Latitude:F8}°, Lon: {Longitude:F8}°";
            }
        }

        /// <summary>
        /// Converts UTM coordinates to Latitude/Longitude (WGS84)
        /// </summary>
        /// <param name="easting">Easting coordinate in meters</param>
        /// <param name="northing">Northing coordinate in meters</param>
        /// <param name="zoneNumber">UTM zone number (1-60)</param>
        /// <param name="isNorthernHemisphere">True for Northern hemisphere, false for Southern</param>
        /// <returns>LatLon object with WGS84 coordinates</returns>
        public static LatLon ConvertUTMToLatLon(double easting, double northing, int zoneNumber, bool isNorthernHemisphere)
        {
            const double k0 = 0.9996; // Scale factor
            double e1 = (1.0 - Math.Sqrt(1.0 - WGS84_E2)) / (1.0 + Math.Sqrt(1.0 - WGS84_E2));

            double x = easting - 500000.0; // Remove false easting
            double y = isNorthernHemisphere ? northing : northing - 10000000.0; // Remove false northing for southern hemisphere

            double M = y / k0;
            double mu = M / (WGS84_A * (1.0 - WGS84_E2 / 4.0 - 3.0 * WGS84_E2 * WGS84_E2 / 64.0 - 5.0 * Math.Pow(WGS84_E2, 3) / 256.0));

            // Footpoint latitude
            double phi1 = mu + (3.0 * e1 / 2.0 - 27.0 * Math.Pow(e1, 3) / 32.0) * Math.Sin(2.0 * mu)
                             + (21.0 * e1 * e1 / 16.0 - 55.0 * Math.Pow(e1, 4) / 32.0) * Math.Sin(4.0 * mu)
                             + (151.0 * Math.Pow(e1, 3) / 96.0) * Math.Sin(6.0 * mu);

            double sinPhi1 = Math.Sin(phi1);
            double cosPhi1 = Math.Cos(phi1);
            double tanPhi1 = Math.Tan(phi1);

            double N1 = WGS84_A / Math.Sqrt(1.0 - WGS84_E2 * sinPhi1 * sinPhi1);
            double T1 = tanPhi1 * tanPhi1;
            double C1 = WGS84_E2 * cosPhi1 * cosPhi1 / (1.0 - WGS84_E2);
            double R1 = WGS84_A * (1.0 - WGS84_E2) / Math.Pow(1.0 - WGS84_E2 * sinPhi1 * sinPhi1, 1.5);
            double D = x / (N1 * k0);

            // Calculate latitude
            double lat = phi1 - (N1 * tanPhi1 / R1) * (D * D / 2.0
                - (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C1 * C1 - 9.0 * WGS84_E2) * Math.Pow(D, 4) / 24.0
                + (61.0 + 90.0 * T1 + 298.0 * C1 + 45.0 * T1 * T1 - 252.0 * WGS84_E2 - 3.0 * C1 * C1) * Math.Pow(D, 6) / 720.0);

            // Calculate longitude
            double lon = (D - (1.0 + 2.0 * T1 + C1) * Math.Pow(D, 3) / 6.0
                + (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C1 * C1 + 8.0 * WGS84_E2 + 24.0 * T1 * T1) * Math.Pow(D, 5) / 120.0) / cosPhi1;

            // Convert to degrees
            lat = lat * 180.0 / Math.PI;
            lon = lon * 180.0 / Math.PI + (zoneNumber - 1) * 6.0 - 180.0 + 3.0;

            return new LatLon { Latitude = lat, Longitude = lon };
        }

        /// <summary>
        /// Converts British National Grid (OSGB36) coordinates to Latitude/Longitude (WGS84)
        /// </summary>
        /// <param name="easting">Easting coordinate in meters</param>
        /// <param name="northing">Northing coordinate in meters</param>
        /// <returns>LatLon object with WGS84 coordinates</returns>
        public static LatLon ConvertBNGToLatLon(double easting, double northing)
        {
            double e2 = 1.0 - (OSGB36_B * OSGB36_B) / (OSGB36_A * OSGB36_A);
            double n = (OSGB36_A - OSGB36_B) / (OSGB36_A + OSGB36_B);

            double lat = OSGB36_LAT0;
            double M = 0.0;

            // Iteratively find the latitude
            while (northing - OSGB36_N0 - M >= 0.00001)
            {
                lat = ((northing - OSGB36_N0 - M) / (OSGB36_A * OSGB36_F0)) + lat;

                double Ma = (1.0 + n + (5.0 / 4.0) * n * n + (5.0 / 4.0) * Math.Pow(n, 3)) * (lat - OSGB36_LAT0);
                double Mb = (3.0 * n + 3.0 * n * n + (21.0 / 8.0) * Math.Pow(n, 3)) * Math.Sin(lat - OSGB36_LAT0) * Math.Cos(lat + OSGB36_LAT0);
                double Mc = ((15.0 / 8.0) * n * n + (15.0 / 8.0) * Math.Pow(n, 3)) * Math.Sin(2.0 * (lat - OSGB36_LAT0)) * Math.Cos(2.0 * (lat + OSGB36_LAT0));
                double Md = (35.0 / 24.0) * Math.Pow(n, 3) * Math.Sin(3.0 * (lat - OSGB36_LAT0)) * Math.Cos(3.0 * (lat + OSGB36_LAT0));

                M = OSGB36_B * OSGB36_F0 * (Ma - Mb + Mc - Md);
            }

            double sinLat = Math.Sin(lat);
            double cosLat = Math.Cos(lat);
            double tanLat = Math.Tan(lat);

            double nu = OSGB36_A * OSGB36_F0 / Math.Sqrt(1.0 - e2 * sinLat * sinLat);
            double rho = OSGB36_A * OSGB36_F0 * (1.0 - e2) / Math.Pow(1.0 - e2 * sinLat * sinLat, 1.5);
            double eta2 = nu / rho - 1.0;

            double VII = tanLat / (2.0 * rho * nu);
            double VIII = tanLat / (24.0 * rho * Math.Pow(nu, 3)) * (5.0 + 3.0 * tanLat * tanLat + eta2 - 9.0 * tanLat * tanLat * eta2);
            double IX = tanLat / (720.0 * rho * Math.Pow(nu, 5)) * (61.0 + 90.0 * tanLat * tanLat + 45.0 * Math.Pow(tanLat, 4));

            double secLat = 1.0 / cosLat;
            double X = secLat / nu;
            double XI = secLat / (6.0 * Math.Pow(nu, 3)) * (nu / rho + 2.0 * tanLat * tanLat);
            double XII = secLat / (120.0 * Math.Pow(nu, 5)) * (5.0 + 28.0 * tanLat * tanLat + 24.0 * Math.Pow(tanLat, 4));

            double dE = easting - OSGB36_E0;
            lat = lat - VII * dE * dE + VIII * Math.Pow(dE, 4) - IX * Math.Pow(dE, 6);
            double lon = OSGB36_LON0 + X * dE - XI * Math.Pow(dE, 3) + XII * Math.Pow(dE, 5);

            // Convert to degrees
            lat = lat * 180.0 / Math.PI;
            lon = lon * 180.0 / Math.PI;

            return new LatLon { Latitude = lat, Longitude = lon };
        }

        static void Main(string[] args)
        {
            Console.WriteLine("=== Easting/Northing to Latitude/Longitude Converter ===\n");

            // Example 1: UTM Conversion
            Console.WriteLine("Example 1: UTM Conversion");
            Console.WriteLine("-------------------------");
            double utmEasting = 651409.903;
            double utmNorthing = 313177.270;
            int zoneNumber = 30;
            bool isNorthern = true;

            LatLon utmResult = ConvertUTMToLatLon(utmEasting, utmNorthing, zoneNumber, isNorthern);
            Console.WriteLine($"Input: UTM Zone {zoneNumber}{(isNorthern ? "N" : "S")}");
            Console.WriteLine($"Easting: {utmEasting:F3} m, Northing: {utmNorthing:F3} m");
            Console.WriteLine($"Output: {utmResult}\n");

            // Example 2: British National Grid Conversion
            Console.WriteLine("Example 2: British National Grid (BNG) Conversion");
            Console.WriteLine("--------------------------------------------------");
            double bngEasting = 530000.0;
            double bngNorthing = 180000.0;

            LatLon bngResult = ConvertBNGToLatLon(bngEasting, bngNorthing);
            Console.WriteLine($"Input: Easting: {bngEasting:F3} m, Northing: {bngNorthing:F3} m");
            Console.WriteLine($"Output: {bngResult}\n");

            // Interactive mode
            Console.WriteLine("\n=== Interactive Mode ===");
            while (true)
            {
                Console.WriteLine("\nSelect coordinate system:");
                Console.WriteLine("1. UTM");
                Console.WriteLine("2. British National Grid (BNG)");
                Console.WriteLine("3. Exit");
                Console.Write("Enter choice (1-3): ");

                string choice = Console.ReadLine();

                if (choice == "3") break;

                try
                {
                    Console.Write("Enter Easting (m): ");
                    double easting = double.Parse(Console.ReadLine());

                    Console.Write("Enter Northing (m): ");
                    double northing = double.Parse(Console.ReadLine());

                    LatLon result;

                    if (choice == "1")
                    {
                        Console.Write("Enter UTM Zone Number (1-60): ");
                        int zone = int.Parse(Console.ReadLine());

                        Console.Write("Northern Hemisphere? (Y/N): ");
                        bool isNorth = Console.ReadLine().ToUpper() == "Y";

                        result = ConvertUTMToLatLon(easting, northing, zone, isNorth);
                    }
                    else if (choice == "2")
                    {
                        result = ConvertBNGToLatLon(easting, northing);
                    }
                    else
                    {
                        Console.WriteLine("Invalid choice!");
                        continue;
                    }

                    Console.WriteLine($"\nResult: {result}");
                    Console.WriteLine($"Google Maps: https://www.google.com/maps?q={result.Latitude:F6},{result.Longitude:F6}");
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"Error: {ex.Message}");
                }
            }

            Console.WriteLine("\nPress any key to exit...");
            Console.ReadKey();
        }
    }
}