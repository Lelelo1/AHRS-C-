using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Spatial.Euclidean;

namespace Logic.Ahrs.Algorithms.Tyrex.Implementations
{
    public class QMichelObs : AttitudeFilter
    {
        double Beta = 0.3;
        double Ka = 2;
        double Km = 1;

        protected override void Update(Vector<double> gyr, Vector<double> acc, Vector<double> mag, double dT)
        {
            updateInternal(gyr, acc, mag, dT);
        }
        void updateInternal(Vector<double> gyr, Vector<double> acc, Vector<double> mag, double dT)
        {
            var q = _Quaternion;

            acc /= Norm(acc);
            mag /= Norm(mag);

            var estimate_A = Quatrotate(q, AccRefNormalized);
            var estimate_M = Quatrotate(q, MagRefNormalized);

            // [acc  mag] matlab list syntax...
            var messure = OneDimensionalMatrix(acc, mag);
            var estimate = OneDimensionalMatrix(estimate_A, estimate_M);

            var delta = 2 * Combined(Ka * Skew(estimate_A), Km * Skew(estimate_M)).Transpose();

            // ((delta * delta' + 1e-5 * eye(3))^-1 * delta)'
            var eye = Matrix<double>.Build.DenseIdentity(3, 3);

            var dq = (messure - estimate) * ((delta * delta.Transpose() + 0.00001 * eye).Inverse() * delta).Transpose();

            var gyrQ = new Quaternion(0, gyr[0], gyr[1], gyr[2]);
            var dqQ = new Quaternion(0, dq[0, 0], dq[0, 1], dq[0, 2]);

            var qDot = 0.5 * (q * gyrQ) + Beta * (q * dqQ);
            q += qDot * dT;
            q /= Norm(q);
            _Quaternion = q;
            Logic.Utils.Log.Message("QMichelObs quaternion became: " + Quaternion);
        }
    }
}
