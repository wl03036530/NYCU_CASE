#include "simulation/kinematics.h"

#include "Eigen/Dense"

#include "acclaim/bone.h"
#include "util/helper.h"
#include <stdio.h>

namespace kinematics {
Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    // TODO
    // You need to return the solution (x) of the linear least squares system:
    //     i.e. find x which min(| Jacobian * x - target |)
    Eigen::Matrix3Xd J = Jacobian.topRows(3);
    Eigen::Vector3d V = target.head<3>();
    Eigen::VectorXd sol = J.transpose() * (J * J.transpose()).inverse() * V;
    //Eigen::VectorXd sol = (Jacobian.transpose() * Jacobian).inverse() * Jacobian.transpose() * target;
    return sol;
    //return Eigen::VectorXd(Jacobian.cols());
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW3 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone, acclaim::Posture& posture)  // target_pos >> ball position
{
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;
    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;
    // TODO
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    size_t bone_num = 1;
    acclaim::Bone*  curr_bone = end_bone;
    while (curr_bone != start_bone && curr_bone != root_bone)
    {
        bone_num += 1;
        curr_bone = curr_bone->parent;
    }
    // HINT:
    // calculate number of bones need to move to perform IK, store in `bone_num`
    // a.k.a. how may bones from end_bone to its parent than to start_bone (include both side)

    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
    Jacobian.setZero();
    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;
        if (desiredVector.norm() < epsilon) {
            break;
        }
        // HINT:
        // Calculate Jacobian, store in `Jacobian`
        curr_bone = end_bone;
        for (size_t i = 0; i < bone_num; i += 1)
        {
            Eigen::Matrix3d rot_mat = curr_bone->rotation.linear();
            Eigen::Vector4d r = curr_bone->start_position;
            Eigen::Vector4d p = end_bone->end_position;
            Eigen::Vector4d delta_position = p - r;
            for (int j = 0; j < 3; j += 1)
            {
                Eigen::Vector3d a = Eigen::Vector3d::Zero();
                if (j == 0)
                    a = (rot_mat* Eigen::Vector3d(1, 0, 0)).normalized();
                else if(j == 1)
                    a = (rot_mat * Eigen::Vector3d(0, 1, 0)).normalized();
                else
                    a = (rot_mat * Eigen::Vector3d(0, 0, 1)).normalized();
                Eigen::Vector3d radius = Eigen::Vector3d::Zero();
                radius << delta_position.head(3);
                Eigen::Vector3d par_dif = a.cross(radius);
                /*printf("%.4lf, %.4lf, %.4lf\n", par_dif.x(), par_dif.y(), par_dif.z());
                system("pause");*/
                Jacobian.col(i * 3 + j) << par_dif, 1;
            }
            curr_bone = curr_bone->parent;
        }
        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);
        // HINT:
        // Change `posture.bone_rotation` based on deltatheta
        curr_bone = end_bone;
        for (size_t i = 0; i < bone_num; i += 1)
        {
            for (int j = 0; j < 3; j += 1)
            {
                posture.bone_rotations[curr_bone->idx][j] += deltatheta[i * 3 + j];
            }
            curr_bone = curr_bone->parent;
        }

    }
    // TODO (Bonus)
    // Return IK is stable?
    // i.e. not swinging its hand in air
    return true;
}
}  // namespace kinematics
