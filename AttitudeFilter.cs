using System;
using System.Threading.Tasks;
using Logic.Ahrs.Service;
using Logic.Services;
using System.Numerics;
using Logic.Models;

// translation of -> https://github.com/tyrex-team/benchmarks-attitude-smartphones/blob/master/src/Filters/AttitudeFilter.m
namespace Logic.Ahrs.Algorithms.Tyrex
{
    // extends 'handle' in matlab source code
    public abstract class AttitudeFilter : IAlgorithm
    {
        public Quaternion Quaternion { get; protected set; } = Quaternion.Identity;
        protected Vector3 AccRef { get; set; } // % Gravity vector in Earth Frame
        protected Vector3 MagRef { get; set; } // % Magnetic vector in Magnetic Earth Frame

        /* not seeming to be used in the AttitudeFilter file or the QMichelObsExtmagWtRep implementation which I looking to use */
        // noises 
        protected double AccMagAngleRef { get; set; } // % Angle between both vectors

        protected double MagRefNorm => MagRef.Length();
        protected Vector3 MagRefNormalized => Vector3.Normalize(MagRef);

        protected double AccRefNorm => AccRef.Length();
        protected Vector3 AccRefNormalized => Vector3.Normalize(AccRef);

        // 'obj' which is passed in n many places seems to be the ase as C#'s 'this'
        protected abstract Quaternion Update(Vector3 gyr, Vector3 acc, Vector3 mag, double dT);

        public bool HasInit { get; set; } = false;
        TaskCompletionSource<Location> HasLocation = new TaskCompletionSource<Location>();
        public void Init(Vector3 gyr, Vector3 acc, Vector3 mag)
        {
            // code found in generateAttitude.m...
            var magCalc = new Geo.Geomagnetism.IgrfGeomagnetismCalculator();// seems to be null always -> new Geo.Geomagnetism.GeomagnetismCalculator();
            
            // location is null when algotithm is started!!!
            var location = new Location(57.7027141, 11.916687); // LocationService.Instance.Location
            var cordinate = new Geo.Coordinate(location.Latitude, location.Longitude);
            var geoMagRes = magCalc.TryCalculate(cordinate, DateTime.UtcNow.Date.AddYears(-7));
            // nanoTesla to microTesla
            var magneticVector = new Vector3((float)geoMagRes.X / 1000, (float)geoMagRes.Y / 1000, (float)geoMagRes.Z / 1000);
            var qTrueToMagnetic = Quaternion.Inverse(Quaternion.CreateFromRotationMatrix(Rotz(geoMagRes.Declination)));
            MagRef = Quatrotate(qTrueToMagnetic, magneticVector);
            var gravityMagnitude = GeoidHeightsDotNet.GeoidHeights.undulation(location.Latitude, location.Latitude); // is -4.5 (is gAtLocation really?)
            var gravityVector = new Vector3(0, 0, (float)-gravityMagnitude);
            AccRef = Quatrotate(qTrueToMagnetic, gravityVector);
            var c = (double)Math.Pow(10, -12); // 1e-12
            MagRef =  MagRef.Abs().Any((double elementValue) => elementValue < c) ? Vector3.Zero : MagRef;
            AccRef = AccRef.Abs().Any((elementValue) => elementValue < c) ? Vector3.Zero : AccRef;

            var H = Vector3.Cross(mag, acc);
            H /= H.Length();

            acc /= acc.Length();
            var M = Vector3.Cross(acc, H);

            // dcm2Quat: The Direction Cosine Matrix to Quaternions block transforms a 3-by-3 direction cosine matrix (DCM) into a four-element unit quaternion vector (q0, q1, q2, q3)
            var R = CreateMatrix(H, M, acc); // enu

            var q = Quaternion.CreateFromRotationMatrix(R);
            Quaternion = q;
            
            HasInit = true;
        }

        public void Update(Reading reading, double dT)
        {
            var r = reading;

            var gyr = new Vector3((float)r.gX, (float)r.gY, (float)r.gZ);
            var acc = new Vector3((float)r.aX, (float)r.aY, (float)r.aZ);
            var mag = new Vector3((float)r.mX, (float)r.mY, (float)r.mZ);

            if (!HasInit)
            {
                Init(gyr, acc, mag);
                Logic.Utils.Log.Message("Algorithm was initialized");
            }

            Update(gyr, acc, mag, dT);
        }

        public Vector3 Quatrotate(Quaternion q, Vector3 v)
        {
            return Vector3.Transform(v, Quaternion.Normalize(q));
        }
        // https://github.com/tyrex-team/benchmarks-attitude-smartphones/blob/master/src/Filters/Implementations/Utils/skew.m
        protected Matrix4x4 Skew(Vector3 v)
        {
            
            // (row, column) https://www.mathworks.com/help/matlab/math/array-indexing.html
            var r12 = -v.Z; // (vector index decending=
            var r13 = v.Y;
            var r23 = -v.X;

            var r21 = v.Z;
            var r31 = -v.Y;
            var r32 = v.X;
            /*
            var arr1 = new double[] { 0, r12, r13 };
            var arr2 = new double[] { r21, 0, r23 };
            var arr3 = new double[] { r31, r32, 0 };
            */
            var m = new Matrix4x4(
                0, r12, r13, 0,
                r21, 0, r23, 0,
                r31, r32, 0, 0,
                0, 0, 0, 1
                );


            return m; 
        }

        protected Matrix4x4 Rotz(double t) // could be wrong!
        {
            t /= (180 * Math.PI);
            var ct = Math.Cos(t);
            var st = Math.Sin(t);

            var v1 = new Vector3((float)ct, (float)-st, (float)0);
            var v2 = new Vector3((float)st, (float)ct, 0);
            var v3 = new Vector3(0, 0, 1);

            return Matrix4x4.CreateWorld(v1, v2, v3);
        }

       
        // practical method https://stackoverflow.com/questions/57921365/why-is-there-no-matrix3x3-in-the-system-numerics-namespace-c-sharp/57923547#comment110976749_57923547
        protected Matrix4x4 CreateMatrix(Vector3 v1, Vector3 v2)
        {
            return new Matrix4x4(
                v1.X, v1.Y, 0, 0,
                v1.Y, v2.Y, 0, 0,
                v1.Z, v2.Z, 0, 0,
                0, 0, 0, 1
                );
        }
        protected Matrix4x4 CreateMatrix(Vector3 v1, Vector3 v2, Vector3 v3)
        {
            return new Matrix4x4(
                v1.X, v2.X, v3.X, 0,
                v1.Y, v2.Y, v3.Y, 0,
                v1.Z, v2.Z, v3.Z, 0,
                0, 0, 0, 1
                );
        }
    }
    enum CordinateSystem
    {
        Uknown,
        ENU,
        NED
    }
    public static class Extension
    {
        public static bool Any(this Vector3 vector, Func<double, bool> condition)
        {
            return condition.Invoke(vector.X) && condition.Invoke(vector.Y) && condition.Invoke(vector.Z);
        }
        public static Vector3 Abs(this Vector3 vector)
        {
            return Vector3.Abs(vector);
        }

        public static Matrix4x4 Add(this Matrix4x4 matrix, double add)
        {
            float a = (float)add;
            matrix.M11 += a;
            matrix.M12 += a;
            matrix.M13 += a;
            matrix.M14 += a;
            matrix.M21 += a;
            matrix.M22 += a;
            matrix.M23 += a;
            matrix.M24 += a;
            matrix.M31 += a;
            matrix.M32 += a;
            matrix.M33 += a;
            matrix.M34 += a;
            matrix.M41 += a;
            matrix.M42 += a;
            matrix.M43 += a;
            matrix.M44 += a;
            return matrix;
        }
        /* conversion between system matrix and mathnet matrix are tested and works */
        /* mathnet matrix.pow could not take -1 exponent
        public static Matrix4x4 Pow(this Matrix4x4 m, int exponent)
        {
            var arr = new float[,] {
                { m.M11, m.M12, m.M13, m.M14 },
                { m.M21, m.M22, m.M23, m.M24 },
                { m.M31, m.M32, m.M33, m.M34 },
                { m.M41, m.M42, m.M43, m.M44 }
            };
            MathNet.Numerics.LinearAlgebra.Matrix<float> mathNetMatrix = MathNet.Numerics.LinearAlgebra.Matrix<float>.Build.DenseOfArray(arr);
            var r = mathNetMatrix.Power(exponent);
            return new Matrix4x4(
                r[0, 0], r[0, 1], r[0, 2], r[0, 3],
                r[1, 0], r[1, 1], r[1, 2], r[1, 3],
                r[2, 0], r[2, 1], r[2, 2], r[2, 3],
                r[3, 0], r[3, 1], r[3, 2], r[3, 3]
                );
        }
        */
        public static Matrix4x4 Inverse(this Matrix4x4 m)
        {
            Matrix4x4 invM;
            Matrix4x4.Invert(m, out invM);
            return invM;
        }
        public static Matrix4x4 Transpose(this Matrix4x4 m)
        {
            return Matrix4x4.Transpose(m);
        }
        public static Quaternion Multiply(this Quaternion q, double value)
        {
            var f = (float)value;
            return Quaternion.Multiply(q, f);
        }
        public static Quaternion Divide(this Quaternion q, double value)
        {
            var f = (float)Math.Pow(value, -1); // as explicit division not being available
            return Quaternion.Multiply(q, f);
        }
        public static Matrix4x4 Multiply(this Matrix4x4 m, double value)
        {
            var f = (float)value;
            return Matrix4x4.Multiply(m, f);
        }
    }
    
}

/* Vector3 and Quaternion length() is same as mathnet matrix L2Norm(), which should
   should be the same as matlab norm() */
