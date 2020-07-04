using System;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Spatial.Euclidean;
// https://github.com/tyrex-team/benchmarks-attitude-smartphones/blob/master/src/Filters/AttitudeFilter.m
namespace Logic.Ahrs.Magnetometer.Tyrex
{
    // extends 'handle' in matlab source code
    public abstract class AttitudeFilter
    {
        public Quaternion Quaternion { get; set; }
        public Vector<double> AccRef { get; set; } // % Gravity vector in Earth Frame
        public Vector<double> MagRef { get; set; } // % Magnetic vector in Magnetic Earth Frame

        /* not seeming to be used in the AttitudeFilter file or the QMichelObsExtmagWtRep implementation which I looking to use */
        // noises 
        protected double AccMagAngleRef { get; set; } // % Angle between both vectors

        protected double MagRefNorm { get; set; }
        protected Vector<double> MagRefNormalized { get; set; }

        protected double AccRefNorm { get; set; }
        protected Vector<double> AccRefNormalized { get; set; }

        // 'obj' which is passed in n many places seems to be the ase as C#'s 'this'
        public abstract Quaternion Update(Vector<double> gyr, Vector<double> acc, Vector<double> mag, double dT);

        public void Init(Vector<double> gyr, Vector<double> acc, Vector<double> mag)
        {
           
            var H = Cross(mag, acc); // 'crossproduct' same as 'cross' ..?
            //H = H / H.Norm(needs arg);
            /* matlab says norm(v) return ecluidian norm https://www.mathworks.com/help/matlab/ref/norm.html
               which should be https://stackoverflow.com/questions/24335521/how-to-get-vector-magnitude-using-mathnet-numerics */
            H /= H.L2Norm(); 

            // set cordinateSystem
            // set quaternion
        }

        public void notifyReferenceVectorChanged()
        {
            MagRefNorm = MagRef.L2Norm();
            AccRefNorm = AccRef.L2Norm();
            AccMagAngleRef = Math.Atan2(Cross(AccRef, MagRef).L2Norm(),
                Vector3D.OfVector(AccRef).DotProduct(Vector3D.OfVector(MagRef)));
            MagRefNormalized = MagRef / MagRefNorm;
            AccRefNormalized = AccRef / AccRefNorm;
        }
        public void GyroInt(Vector<double> gyr, Vector<double> acc, double dT)
        {
            var q = Quaternion;

            q = q * new Quaternion(0.5 * dT * DenseVector.OfVector(gyr));
            q /= q.Norm;

            Quaternion = q;
        }
        // helper method to give `cross` to LinearAlgebra namespace: https://stackoverflow.com/questions/11759720/cross-product-using-math-net-numerics-with-c-sharp
        public Vector<double> Cross(Vector<double> left, Vector<double> right)
        {
            if ((left.Count != 3 || right.Count != 3))
            {
                string message = "Vectors must have a length of 3. Tyrex.cs";
                throw new Exception(message);
            }
            var result = new DenseVector(3);

            result[0] = left[1] * right[2] - left[2] * right[1];
            result[1] = -left[0] * right[2] + left[2] * right[0];
            result[2] = left[0] * right[1] - left[1] * right[0];

            return result;
        }

        public Vector<double> Quatrotate(Quaternion q, Vector<double> v)
        {
            q = q.Normalized;

            var qX = (float)q.ImagX;
            var qY = (float)q.ImagY;
            var qZ = (float)q.ImagZ;
            var qW = (float)q.Real;
            var systemQ = new System.Numerics.Quaternion(new System.Numerics.Vector3(qX, qY, qZ), qW);

            var vX = (float)v[0];
            var vY = (float)v[1];
            var vZ = (float)v[2];
            var systemV = new System.Numerics.Vector3(vX, vY, vZ);

            var r = System.Numerics.Vector3.Transform(systemV, systemQ);

            return Vector<double>.Build.DenseOfArray(new double[] { r.X, r.Y, r.Z });
        }
    }
    enum CordinateSystem
    {
        Uknown,
        ENU,
        NED
    }

    
}
