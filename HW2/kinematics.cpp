#include "simulation/kinematics.h"

#include "Eigen/Dense"

#include "acclaim/bone.h"
#include "util/helper.h"

void BFS(acclaim::Bone* bone, bool visited[], const acclaim::Posture& posture)
{
    visited[bone->idx] = true;
    bone->start_position = bone->parent->end_position;

    bone->rotation = bone->parent->rotation * bone->rot_parent_current;
    bone->rotation *= util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);

    bone->end_position = bone->start_position;
    bone->end_position += bone->rotation * (bone->dir.normalized() * bone->length);

    acclaim::Bone* temp = bone->sibling;
    while (temp != nullptr)
    {
        if ( !visited[temp->idx])
        {
            BFS(temp, visited, posture);
        }
        temp = temp->sibling;
    }
    if (bone->child != nullptr)
        if (!visited[bone->child->idx])
            BFS(bone->child, visited, posture);
}

namespace kinematics {
void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO
    // This function will be called with bone == root bone of the skeleton
    // You should set these variables:
    //     bone->start_position = Eigen::Vector4d::Zero();
    //     bone->end_position = Eigen::Vector4d::Zero();
    //     bone->rotation = Eigen::Matrix4d::Zero();
    // The sample above just set everything to zero

    bool visited[30] = {};
    // root
    bone->start_position = posture.bone_translations[0];

    bone->rotation = bone->rot_parent_current;
    bone->rotation = util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);

    bone->end_position = bone->start_position;
    bone->end_position += bone->rotation * bone->dir * bone->length;
    //printf("%.3lf, %.3lf, %.3lf, %.3lf\n", posture.bone_rotations[bone->idx].x(), posture.bone_rotations[bone->idx].y(), posture.bone_rotations[bone->idx].z(), posture.bone_rotations[bone->idx].w());
    visited[0] = true;

    BFS(bone->child, visited, posture);
}

std::vector<acclaim::Posture> timeWarper(const std::vector<acclaim::Posture>& postures, int keyframe_old,  int keyframe_new)    // 160 150
{
    int total_frames = static_cast<int>(postures.size());
    int total_bones = static_cast<int>(postures[0].bone_rotations.size());
    std::vector<acclaim::Posture> new_postures = postures;

    for (int i = 0; i < total_frames; ++i) {
        for (int j = 0; j < total_bones; ++j) {
            // TODO
            // You should set these variables:
            //     new_postures[i].bone_translations[j] = postures[i].bone_translations[j];
            //     new_postures[i].bone_rotations[j] = postures[i].bone_rotations[j];
            // The sample above just change nothing

            if (i <= keyframe_new)
            {
                double interval = keyframe_old*1.0 / keyframe_new*1.0;
                double now_idx = interval * i;
                
                int Floor = now_idx;
                int Ceiling = now_idx + 1;

                double dist_to_p1 = now_idx - now_idx;

                Eigen::Quaterniond transform;

                transform.x() = postures[Floor].bone_translations[j].x();
                transform.y() = postures[Floor].bone_translations[j].y();
                transform.z() = postures[Floor].bone_translations[j].z();
                transform.w() = postures[Floor].bone_translations[j].w();

                Eigen::Quaterniond transform_next;
                transform_next.x() = postures[Ceiling].bone_translations[j].x();
                transform_next.y() = postures[Ceiling].bone_translations[j].y();
                transform_next.z() = postures[Ceiling].bone_translations[j].z();
                transform_next.w() = postures[Ceiling].bone_translations[j].w();

                Eigen::Quaterniond rotation;
                rotation.x() = postures[Floor].bone_rotations[j].x();
                rotation.y() = postures[Floor].bone_rotations[j].y();
                rotation.z() = postures[Floor].bone_rotations[j].z();
                rotation.w() = postures[Floor].bone_rotations[j].w();

                Eigen::Quaterniond rotation_next;
                rotation_next.x() = postures[Ceiling].bone_rotations[j].x();
                rotation_next.y() = postures[Ceiling].bone_rotations[j].y();
                rotation_next.z() = postures[Ceiling].bone_rotations[j].z();
                rotation_next.w() = postures[Ceiling].bone_rotations[j].w();

                Eigen::Quaterniond new_transform = transform.slerp(dist_to_p1, transform_next);
                Eigen::Quaterniond new_rotation = rotation.slerp(dist_to_p1, rotation_next);

                new_postures[i].bone_translations[j] = Eigen::Vector4d(new_transform.x(), new_transform.y(), new_transform.z(), new_transform.w());
                new_postures[i].bone_rotations[j] = Eigen::Vector4d(new_rotation.x(), new_rotation.y(), new_rotation.z(), new_rotation.w());
            }
            else
            {
                double interval = (total_frames - keyframe_old) * 1.0 / (total_frames - keyframe_new) * 1.0;
                double now_idx = interval * (i - keyframe_new) + keyframe_old -1;
                int Floor = now_idx;
                int Ceiling = now_idx + 1;

                double dist_to_p1 = now_idx - Floor;

                Eigen::Quaterniond transform;

                transform.x() = postures[Floor].bone_translations[j].x();
                transform.y() = postures[Floor].bone_translations[j].y();
                transform.z() = postures[Floor].bone_translations[j].z();
                transform.w() = postures[Floor].bone_translations[j].w();

                Eigen::Quaterniond transform_next;
                transform_next.x() = postures[Ceiling].bone_translations[j].x();
                transform_next.y() = postures[Ceiling].bone_translations[j].y();
                transform_next.z() = postures[Ceiling].bone_translations[j].z();
                transform_next.w() = postures[Ceiling].bone_translations[j].w();

                Eigen::Quaterniond rotation;
                rotation.x() = postures[Floor].bone_rotations[j].x();
                rotation.y() = postures[Floor].bone_rotations[j].y();
                rotation.z() = postures[Floor].bone_rotations[j].z();
                rotation.w() = postures[Floor].bone_rotations[j].w();

                Eigen::Quaterniond rotation_next;
                rotation_next.x() = postures[Ceiling].bone_rotations[j].x();
                rotation_next.y() = postures[Ceiling].bone_rotations[j].y();
                rotation_next.z() = postures[Ceiling].bone_rotations[j].z();
                rotation_next.w() = postures[Ceiling].bone_rotations[j].w();

                Eigen::Quaterniond new_transform = transform.slerp(dist_to_p1, transform_next);
                Eigen::Quaterniond new_rotation = rotation.slerp(dist_to_p1, rotation_next);

                new_postures[i].bone_translations[j] = Eigen::Vector4d(new_transform.x(), new_transform.y(), new_transform.z(), new_transform.w());
                new_postures[i].bone_rotations[j] = Eigen::Vector4d(new_rotation.x(), new_rotation.y(), new_rotation.z(), new_rotation.w());
            }
        }
    }
    return new_postures;
}
}  // namespace kinematics
