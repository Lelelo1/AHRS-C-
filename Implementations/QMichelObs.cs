using System;
using System.Numerics;

namespace Logic.Ahrs.Algorithms.Tyrex.Implementations
{
    public class QMichelObs : AttitudeFilter
    {
        double Beta = 0.3;
        double Ka = 2;
        double Km = 1;

        protected override Quaternion Update(Vector3 gyr, Vector3 acc, Vector3 mag, double dT)
        {
            updateInternal(gyr, acc, mag, dT);
            return Quaternion;
        }
        void updateInternal(Vector3 gyr, Vector3 acc, Vector3 mag, double dT)
        {
            var q = Quaternion;

            acc /= acc.Length();
            mag /= mag.Length();

            var estimate_A = Quatrotate(q, AccRefNormalized);
            var estimate_M = Quatrotate(q, MagRefNormalized);

            var messure = CreateMatrix(acc, mag);
            var estimate = CreateMatrix(estimate_A, estimate_M);

            var sA = Skew(estimate_A).Multiply(Ka);
            var sM = Skew(estimate_M).Multiply(Km);
            var trans = Matrix4x4.Add(sA, sM).Transpose();
            var delta = trans.Multiply(2);

            // ((delta * delta' + 1e-5 * eye(3))^-1 * delta)'
            var eye = Matrix4x4.Identity.Add(0.00001);
            var deltas = ((delta * delta.Transpose() * eye).Inverse() * delta).Transpose();
            var dq = (messure - estimate) * deltas;

            var gyrQ = new Quaternion(gyr.X, gyr.Y, gyr.Z, 0);
            
            var dqQ = Quaternion.CreateFromRotationMatrix(dq);
            dqQ.W = 0;
            
            //var dqQ = new Quaternion(dq.M11, dq.M21, dq.M31, 0);

            var qDot = (q * gyrQ).Multiply(0.5) + (q * dqQ).Multiply(Beta);
            q += qDot.Multiply(dT);
            q = q.Divide(q.Length());

            Quaternion = q;
            Logic.Utils.Log.Message("QMichelObs quaternion became: " + Quaternion);
        }
    }
}
