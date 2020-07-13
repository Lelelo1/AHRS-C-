using System;
using System.Threading.Tasks;
using Logic.Ahrs.Service;
using Logic.Models;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;
using MathNet.Spatial.Units;

using System.Linq;

// translation of -> https://github.com/tyrex-team/benchmarks-attitude-smartphones/blob/master/src/Filters/AttitudeFilter.m
namespace Logic.Ahrs.Algorithms.Tyrex
{
    // extends 'handle' in matlab source code
    public abstract class AttitudeFilter : IAlgorithm
    {
        public System.Numerics.Quaternion Quaternion => (ToMagnetic * _Quaternion).AsSystemQuaternion();

        protected Quaternion _Quaternion { get; set; }
        protected Vector<double> AccRef { get; set; } // % Gravity vector in Earth Frame
        protected Vector<double> MagRef { get; set; } // % Magnetic vector in Magnetic Earth Frame

        /* not seeming to be used in the AttitudeFilter file or the QMichelObsExtmagWtRep implementation which I looking to use */
        // noises 
        protected double AccMagAngleRef { get; set; } // % Angle between both vectors

        protected double MagRefNorm { get; set; }
        protected Vector<double> MagRefNormalized { get; set; }

        protected double AccRefNorm { get; set; }
        protected Vector<double> AccRefNormalized { get; set; }


        protected abstract void Update(Vector<double> gyr, Vector<double> acc, Vector<double> mag, double dT);

        public bool HasInit { get; set; } = false;
        TaskCompletionSource<Location> HasLocation = new TaskCompletionSource<Location>();

        Quaternion ToTrue { get; set; }
        Quaternion ToMagnetic { get; set; }

        protected CordinateSystem CordinateSystem { get; set; } = CordinateSystem.ENU;

        
        void Setup()
        {

            Console.WriteLine("doing setup");

            // createContextFromLocationAndDate.m
            var location = new Location(57.7027141, 11.916687); // LocationService.Instance.Location
            var date = DateTime.UtcNow;//DateTime.UtcNow.Date.AddYears(-7);

            var toNED = DCMToQuat(RotY(180) * RotZ(90));
            var toENU = QuatInv(toNED);

            var magCalc = new Geo.Geomagnetism.WmmGeomagnetismCalculator();
            var cordinate = new Geo.Coordinate(location.Latitude, location.Longitude);
            var geoMagRes = magCalc.TryCalculate(cordinate, date);

            // %Magnetic field context & Transform nanoTesla to microTesla
            var mX = geoMagRes.X / 1000;
            var mY = geoMagRes.Y / 1000;
            var mZ = geoMagRes.Z / 1000;

            var magneticVector = new double[] { mX, mY, mZ }.ToVector();
            // mine %v = [-3.32527858581955, 16.156762327648, -48.2670628388773];
            // test matlab 0.8465   22.8248  -41.3634
            var magneticDeclination = geoMagRes.Declination;

            if (CordinateSystem == CordinateSystem.ENU)
            {
                magneticVector = Quatrotate(toENU, magneticVector);
                magneticDeclination = -magneticDeclination;
            }
            
            // % Gravity context
            var gravityMagnitude = 9.81; // don't differ much and can be set arbitrarily in most use cases, can otherwise be set based on location
            var gravityVector = new double[] { 0, 0, -gravityMagnitude }.ToVector();

            if(CordinateSystem == CordinateSystem.ENU)
            {
                gravityVector = Quatrotate(toENU, gravityVector);
            }

            //magnetic declination: Mine -4.19322488655381. Matlab -2.1241

            magneticVector = new double[] { 0.8465, 22.8248, -41.3634 }.ToVector();
            // start generateAttitude.m
            ToTrue = DCMToQuat(Rotz(magneticDeclination)); // 0.999330236733822 0 0 0.036593556739323
            ToMagnetic = QuatInv(ToTrue); // 0.999 0 0 -0.0366122580820162

            // MagRef = new double[] { -0.000, 22.8405, -41.3634 }.ToVector(); // more or less the same
            MagRef = Quatrotate(ToMagnetic, magneticVector);
            AccRef = Quatrotate(ToMagnetic, gravityVector);
            // MagRef[0] = 0; // almost same but creates a drift leftwards

            Func<double, bool> where = (double value) => Math.Abs(value) < Math.Pow(10, -12);
            MagRef = MagRef.SetAny(where, 0);
            AccRef = AccRef.SetAny(where, 0);

            // MagRef = [-0.8239503291328, 22.8256251424549, -41.3634]

            var mValue = MagRef.BoolVector((value) => Math.Abs(value) > 0).Sum();
            var aValue = AccRef.BoolVector((value) => Math.Abs(value) > 0).Sum();
            if(mValue != 2 || aValue != 1)
            {
                Console.WriteLine("Warning, Reference vectors are not well constructed, AttitudeFilter.cs"); // % This should never happen
            }
            
            ReferenceVectorChanged();
        }

        public void Init(Vector<double> gyr, Vector<double> acc, Vector<double> mag)
        {
            Setup();
            // end generateAttitude.m

            var H = Cross(mag, acc);
            H /= Norm(H);

            acc /= Norm(acc);
            var M = Cross(acc, H);

            acc = CordinateSystem == CordinateSystem.NED ? -acc : acc;
            var R = CreateMatrix(H, M, acc, vectorsAsColumns: true); 

            var q = DCMToQuat(R);
            _Quaternion = q;
            
            HasInit = true;
        }

        public void Update(Reading reading, double dT)
        {
            var r = reading;

            var gyr = new double[] { r.gX, r.gY, r.gZ }.ToVector();
            var acc = new double[] { r.aX, r.aY, r.aZ }.ToVector();
            var mag = new double[] { r.mX, r.mY, r.mZ }.ToVector();

            if (!HasInit)
            {
                Init(gyr, acc, mag);
                Logic.Utils.Log.Message("Algorithm was initialized");
            }
           
            Update(gyr, acc, mag, dT);
        }

        void ReferenceVectorChanged()
        {
            MagRefNorm = Norm(MagRef);
            AccRefNorm = Norm(AccRef);
            AccMagAngleRef = Math.Atan2(Norm(Cross(AccRef, MagRef)), Dot(AccRef, MagRef)); // not used anywhere currently

            MagRefNormalized = MagRef / MagRefNorm;
            AccRefNormalized = AccRef / AccRefNorm;
        }

        // is ok!
        protected Quaternion QuatInv(Quaternion q)
        {
            return new Quaternion(q.Real, -q.ImagX, -q.ImagY, -q.ImagZ);
        }

        // is ok!
        protected Matrix<double> Combined(Matrix<double> m1, Matrix<double> m2)
        {
            /* [m1]
               [m2] */
            // https://stackoverflow.com/questions/50902794/matrix-concatenation-using-mathnet-numerics-linearalgebra-in-net-framework/62815400#62815400
            return Matrix<double>.Build.DenseOfMatrixArray(new Matrix<double>[,] { { m1 }, { m2 } });
        }

        // is ok!
        protected static double Dot(Vector<double> v1, Vector<double> v2)
        {
            return v1.DotProduct(v2);
        }

        // ok!
        protected Vector<double> Quatrotate(Quaternion q, Vector<double> v, string form = "short")
        {
            var qw = q.Real; var qx = q.ImagX; var qy = q.ImagY; var qz = q.ImagZ;
            var vx = v[0]; var vy = v[1]; var vz = v[2];

            Vector<double> vo = null;

            if (form == "long")
            {
                // seperate expressions for easier debugging
                var xFirst = vx * (qw.P(2) + qx.P(2) - qy.P(2) - qz.P(2));
                var xSecond = vy * (2 * qw * qz + 2 * qx * qy);
                var xThird = vz * (2 * qw * qy - 2 * qx * qz);
                var x = xFirst + xSecond - xThird;

                var yFirst = vy * (qw.P(2) - qx.P(2) + qy.P(2) - qz.P(2));
                var ySecond = vx * (2 * qw * qz - 2 * qx * qy);
                var yThird = vz * (2 * qw * qx + 2 * qy * qz);
                var y = yFirst - ySecond  + yThird;

                var zFirst = vz * (qw.P(2) - qx.P(2) - qy.P(2) + qz.P(2));
                var zSecond = vx * (2 * qw * qy + 2 * qx * qz);
                var zThird = vy * (2 * qw * qx - 2 * qy * qz);
                var z = zFirst + zSecond - zThird;

                vo = new double[] { x, y, z }.ToVector();
            }
            else if (form == "short")
            {
                var x = vy * (2 * qw * qz + 2 * qx * qy) - vx * (2 * qy.P(2) + 2 * qz.P(2) - 1) - vz * (2 * qw * qy - 2 * qx * qz);
                var y = vz * (2 * qw * qx + 2 * qy * qz) - vx * (2 * qw * qz - 2 * qx * qy) - vy * (2 * qx.P(2) + 2 * qz.P(2) - 1);
                var z = vx * (2 * qw * qy + 2 * qx * qz) - vz * (2 * qx.P(2) + 2 * qy.P(2) - 1) - vy * (2 * qw * qx - 2 * qy * qz);

                vo = new double[] { x, y, z }.ToVector();
            }

            return vo;
        }

        protected Matrix<double> RotY(double degs)
        {
            return Matrix3D.RotationAroundYAxis(Angle.FromDegrees(degs));
        }

        protected Matrix<double> RotZ(double degs)
        {
            return Matrix3D.RotationAroundZAxis(Angle.FromDegrees(degs));
        }

        // ok!
        protected double Norm(Vector<double> v)
        {
            return v.L2Norm();
        }

        // ok!
        protected double Norm(Quaternion q)
        {
            return q.Norm;
        }

        // ok!
        protected Vector<double> Cross(Vector<double> v1, Vector<double> v2)
        {
            return Vector3D.OfVector(v1).CrossProduct(Vector3D.OfVector(v2)).ToVector();
        }

        protected Quaternion DCMToQuat(Matrix<double> dcm)
        {
            var tr = dcm[0, 0] + dcm[1, 1] + dcm[2, 2];
            Quaternion quat = MathNet.Spatial.Euclidean.Quaternion.Zero;

            bool cond1 = tr > 0; // ok!
            bool cond2 = (dcm[0, 0] > dcm[1, 1]) && (dcm[0, 0] > dcm[2, 2]); // ok!
            bool cond3 = (dcm[1, 1] > dcm[2, 2]); // ok!
            // else ok!

            if (cond1)
            {
                var s = Math.Sqrt(tr + 1) * 2;

                quat = new Quaternion(
                    0.25 * s,
                    (dcm[1, 2] - dcm[2, 1]) / s,
                    (dcm[2, 0] - dcm[0, 2]) / s,
                    (dcm[0, 1] - dcm[1, 0]) / s
                    );
            }
            else if (cond2)
            {
                var s = Math.Sqrt(1 + dcm[0, 0] - dcm[1, 1] - dcm[2, 2]) * 2;

                quat = new Quaternion(
                    (dcm[1, 2] - dcm[2, 1]) / s,
                    0.25 * s,
                    (dcm[1, 0] + dcm[0, 1]) / s,
                    (dcm[2, 0] + dcm[0, 2]) / s
                    );
            }
            else if (cond3)
            {
                var s = Math.Sqrt(1 + dcm[1, 1] - dcm[0, 0] - dcm[2, 2]) * 2;
                quat = new Quaternion(
                    (dcm[2, 0] - dcm[0, 2]) / s,
                    (dcm[1, 0] + dcm[0, 1]) / s,
                    0.25 * s,
                    (dcm[2, 1] + dcm[1, 2]) / s
                    );
            }
            else
            {
                var s = Math.Sqrt(1 + dcm[2, 2] - dcm[0, 0] - dcm[1, 1]) * 2;

                quat = new Quaternion(
                    (dcm[0, 1] - dcm[1, 0]) / s,
                    (dcm[2, 0] + dcm[0, 2]) / s,
                    (dcm[2, 1] + dcm[1, 2]) / s,
                    0.25 * s
                    );
            }
            return quat;
        }

        // systen numerics transform and quatrotate.m give different result, can't test if tyrex translated quatrotate is correct
        /*
        public Vector3 MatlabQuatrotate(Quaternion q, Vector3 v)
        {
            // case short
            // P is power
            var qw = q.W; var qx = q.X; var qy = q.Y; var qz = q.Z;
            var vx = v.X; var vy = v.Y; var vz = v.Z;

            var x = vy*(2*qw*qz + 2*qx*qy) - vx*(2*qy.P(2) + 2*qz.P(2)-1) - vz*(2*qw*qy - 2*qx*qz);
            var y = vz*(2*qw*qx + 2*qy*qz) - vx*(2*qw*qz - 2*qx*qy) - vy*(2*qx.P(2) + 2*qz.P(2) - 1);
            var z = vx*(2*qw*qy + 2*qx*qz) - vz*(2*qx.P(2) + 2*qy.P(2) - 1) - vy*(2*qw*qx - 2*qy*qz);

            return new Vector3(x, y, z);
        }
        */

        // skew printed with skew.m with vector [2, 4, 8] 
        //  which gives following matrix:       0    -8     4
        //                                      8     0    -2
        //                                     -4     2     0
        protected Matrix<double> Skew(Vector<double> v)
        {
            // (row, column) https://www.mathworks.com/help/matlab/math/array-indexing.html
            var r12 = -v[2]; // vector index decending
            var r13 = v[1];
            var r23 = -v[0];

            var r21 = v[2];
            var r31 = -v[1];
            var r32 = v[0];

            var arr1 = new double[] { 0, r12, r13 };
            var arr2 = new double[] { r21, 0, r23 };
            var arr3 = new double[] { r31, r32, 0 };

            return CreateMatrix(arr1.ToVector(), arr2.ToVector(), arr3.ToVector());
            // Tested is correct!
        }

        protected Matrix<double> Rotz(double t) // could be wrong!
        {
            t = (t / 180) * Math.PI;

            var ct = (float)Math.Cos(t);
            var st = (float)Math.Sin(t);

            var arr1 = new double[] { ct, -st, 0 };
            var arr2 = new double[] { st, ct, 0 };
            var arr3 = new double[] { 0, 0, 1 };

            return CreateMatrix(arr1.ToVector(), arr2.ToVector(), arr3.ToVector());
            // Tested is correct!
        }

        // placing vectors in rows, which flows better with the matlab syntax: [ r11, r12, r13; 
        //                                                                       r21, r22, r23 ]
        // and is more left to right readable
        protected Matrix<double> CreateMatrix(Vector<double> v1, Vector<double> v2)
        {
            return Matrix<double>.Build.DenseOfRowVectors(v1, v2);
        }
        protected Matrix<double> CreateMatrix(Vector<double> v1, Vector<double> v2, Vector<double> v3, bool vectorsAsColumns = false)
        {
            return vectorsAsColumns ? Matrix<double>.Build.DenseOfColumnVectors(v1, v2, v3) : Matrix<double>.Build.DenseOfRowVectors(v1, v2, v3);
        }

        protected Matrix<double> OneDimensionalMatrix(Vector<double> v1, Vector<double> v2)
        {
            var arr1 = v1.AsArray();
            var arr2 = v2.AsArray();
            return Matrix<double>.Build.DenseOfRowArrays(new double[] {
                arr1[0], arr1[1], arr1[2],
                arr2[1], arr2[2],arr2[2]
            });
        }
    }
    public enum CordinateSystem
    {
        ENU,
        NED
    }
    public static class Extension
    {
        public static Vector<double> ToVector(this double[] array)
        {
            return Vector<double>.Build.Dense(array);
        }
        public static float P(this float f, int p)
        {
            return (float)Math.Pow(f, p);
        }
        public static double P(this double d, int p)
        {
            return Math.Pow(d, p);
        }


        public static Vector<double> SetAny(this Vector<double> v, Func<double, bool> condition, double value)
        {

            Func<int, double> setter = (i) =>
            {
                var any = v[i];
                return condition.Invoke(any) ? value : any;
            };

            return SetVector(v, setter);
        } // abs(v) > 0 matlab syntax, returns 1 or 0 on each element
        public static Vector<double> BoolVector(this Vector<double> v, Func<double, bool> condition)
        {
            var newV = new double[] { v[0], v[1], v[2] }.ToVector();

            Func<int, double> setter = (i) =>
            {
                var any = newV[i];
                return condition.Invoke(any) ? 1 : 0;
            };
            return SetVector(newV, setter);
        }
        static Vector<double> SetVector(Vector<double> v, Func<int, double> set)
        {
            for (int i = 0; i < v.Count; i++)
            {
                v[i] = set.Invoke(i);
            };
            return v;
        }


        public static System.Numerics.Quaternion AsSystemQuaternion(this Quaternion q)
        {
            return new System.Numerics.Quaternion((float)q.ImagX, (float)q.ImagY, (float)q.ImagZ, (float)q.Real);
        }
    }
    
}



// I also noticed that `filter.MagRef = quatrotate(qTrueToMagnetic, context.magnetic.vector);` operation creates a `x = 0` in the vector which is necessary to pass the check()