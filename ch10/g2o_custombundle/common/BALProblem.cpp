#include "BALProblem.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tools/random.h"
#include "tools/rotation.h"

typedef Eigen::Map<Eigen::VectorXd> VectorRef;
typedef Eigen::Map<const Eigen::VectorXd> ConstVectorRef;

template <typename T>
void FscanfOrDie(FILE *fptr, const char *format, T *value)
{
    int num_scanned = fscanf(fptr, format, value);
    if (num_scanned != 1)
        std::cerr << "Invalid UW data file. ";
}

void PerturbPoint3(const double sigma, double *point)
{
    for (int i = 0; i < 3; ++i)
        point[i] += RandNormal() * sigma;
}

double Median(std::vector<double> *data)
{
    int n = data->size();
    std::vector<double>::iterator mid_point = data->begin() + n / 2;
    std::nth_element(data->begin(), mid_point, data->end());
    return *mid_point;
}

BALProblem::BALProblem(const std::string &filename, bool use_quaternions)
{
    FILE *fptr = fopen(filename.c_str(), "r");

    if (fptr == NULL)
    {
        std::cerr << "Error: unable to open file " << filename;
        return;
    };

    // This wil die horribly on invalid files. Them's the breaks.
    // 相机数量(16)，路标数量(22106)，总的观测数量(83718)
    FscanfOrDie(fptr, "%d", &num_cameras_);
    FscanfOrDie(fptr, "%d", &num_points_);
    FscanfOrDie(fptr, "%d", &num_observations_);

    std::cout << "Header: " << num_cameras_
              << " " << num_points_
              << " " << num_observations_;

    // 分配空间 sizeof(int) * 83718
    point_index_ = new int[num_observations_];
    // 分配空间 sizeof(int) * 83718
    camera_index_ = new int[num_observations_];
    // 分配空间 sizeof(int) * 2 * 83718
    // 第0个位置存储相机id，第1个位置存储路标点id , (2 , 3)存储观察到的像素点坐标
    observations_ = new double[2 * num_observations_];

    // 分配空间 9 * 16 + 3 * 22106
    // (0 ,143) 存储16个相机参数信息，每个相机参数信息包含9个参数，分别表示为3个旋转向量，3个平移向量，焦距1个，2个畸变系数
    // （144, ～） 存储路标点的在世界坐标系下的三维点
    num_parameters_ = 9 * num_cameras_ + 3 * num_points_;
    parameters_ = new double[num_parameters_];

    for (int i = 0; i < num_observations_; ++i)
    {
        FscanfOrDie(fptr, "%d", camera_index_ + i);
        FscanfOrDie(fptr, "%d", point_index_ + i);
        for (int j = 0; j < 2; ++j)
        {
            FscanfOrDie(fptr, "%lf", observations_ + 2 * i + j);
        }
    }

    for (int i = 0; i < num_parameters_; ++i)
    {
        FscanfOrDie(fptr, "%lf", parameters_ + i);
    }

    fclose(fptr);

    use_quaternions_ = use_quaternions;
    if (use_quaternions)
    {
        // Switch the angle-axis rotations to quaternions.
        // 四元数 4 + 平移 3 + 相机内参 3 = 10
        num_parameters_ = 10 * num_cameras_ + 3 * num_points_;
        double *quaternion_parameters = new double[num_parameters_];
        double *original_cursor = parameters_;
        double *quaternion_cursor = quaternion_parameters;
        for (int i = 0; i < num_cameras_; ++i)
        {
            AngleAxisToQuaternion(original_cursor, quaternion_cursor);
            quaternion_cursor += 4;
            original_cursor += 3;
            for (int j = 4; j < 10; ++j)
            {
                *quaternion_cursor++ = *original_cursor++;
            }
        }
        // Copy the rest of the points.
        for (int i = 0; i < 3 * num_points_; ++i)
        {
            *quaternion_cursor++ = *original_cursor++;
        }
        // Swap in the quaternion parameters.
        delete[] parameters_;
        parameters_ = quaternion_parameters;
    }
}

void BALProblem::WriteToFile(const std::string &filename) const
{
    FILE *fptr = fopen(filename.c_str(), "w");

    if (fptr == NULL)
    {
        std::cerr << "Error: unable to open file " << filename;
        return;
    }

    fprintf(fptr, "%d %d %d %d\n", num_cameras_, num_cameras_, num_points_, num_observations_);

    for (int i = 0; i < num_observations_; ++i)
    {
        fprintf(fptr, "%d %d", camera_index_[i], point_index_[i]);
        for (int j = 0; j < 2; ++j)
        {
            fprintf(fptr, " %g", observations_[2 * i + j]);
        }
        fprintf(fptr, "\n");
    }

    for (int i = 0; i < num_cameras(); ++i)
    {
        double angleaxis[9];
        if (use_quaternions_)
        {
            //OutPut in angle-axis format.
            QuaternionToAngleAxis(parameters_ + 10 * i, angleaxis);
            memcpy(angleaxis + 3, parameters_ + 10 * i + 4, 6 * sizeof(double));
        }
        else
        {
            memcpy(angleaxis, parameters_ + 9 * i, 9 * sizeof(double));
        }
        for (int j = 0; j < 9; ++j)
        {
            fprintf(fptr, "%.16g\n", angleaxis[j]);
        }
    }

    const double *points = parameters_ + camera_block_size() * num_cameras_;
    for (int i = 0; i < num_points(); ++i)
    {
        const double *point = points + i * point_block_size();
        for (int j = 0; j < point_block_size(); ++j)
        {
            fprintf(fptr, "%.16g\n", point[j]);
        }
    }

    fclose(fptr);
}

// Write the problem to a PLY file for inspection in Meshlab or CloudCompare
void BALProblem::WriteToPLYFile(const std::string &filename) const
{
    std::ofstream of(filename.c_str());

    of << "ply"
       << '\n'
       << "format ascii 1.0"
       << '\n'
       << "element vertex " << num_cameras_ + num_points_
       << '\n'
       << "property float x"
       << '\n'
       << "property float y"
       << '\n'
       << "property float z"
       << '\n'
       << "property uchar red"
       << '\n'
       << "property uchar green"
       << '\n'
       << "property uchar blue"
       << '\n'
       << "end_header" << std::endl;

    // Export extrinsic data (i.e. camera centers) as green points.
    double angle_axis[3];
    double center[3];
    for (int i = 0; i < num_cameras(); ++i)
    {
        const double *camera = cameras() + camera_block_size() * i;
        CameraToAngelAxisAndCenter(camera, angle_axis, center);
        std::cout << center[0] << " , " << center[1] << " , " << center[2] << std::endl;
        of << center[0] << ' ' << center[1] << ' ' << center[2]
           << " 0 255 0" << '\n';
    }

    // Export the structure (i.e. 3D Points) as white points.
    const double *points = parameters_ + camera_block_size() * num_cameras_;
    for (int i = 0; i < num_points(); ++i)
    {
        const double *point = points + i * point_block_size();
        for (int j = 0; j < point_block_size(); ++j)
        {
            of << point[j] << ' ';
        }
        of << "255 255 255\n";
    }
    of.close();
}

void BALProblem::CameraToAngelAxisAndCenter(const double *camera,
                                            double *angle_axis,
                                            double *center) const
{
    VectorRef angle_axis_ref(angle_axis, 3);
    if (use_quaternions_)
    {
        QuaternionToAngleAxis(camera, angle_axis);
    }
    else
    {
        angle_axis_ref = ConstVectorRef(camera, 3);
    }

    // c = -R't
    // Pc = RPw + T
    // Pw = R'(Pc - T)
    // 已知Pc = (0,0,0)T
    // Pw = -R'T
    // R表示世界坐标系到相机坐标系的旋转矩阵
    // T表示平移量
    // 轴角取负号表示矩阵的逆
    Eigen::VectorXd inverse_rotation = -angle_axis_ref;
    AngleAxisRotatePoint(inverse_rotation.data(),
                         camera + camera_block_size() - 6,
                         center);

    VectorRef(center, 3) *= -1.0;
}

void BALProblem::AngleAxisAndCenterToCamera(const double *angle_axis,
                                            const double *center,
                                            double *camera) const
{
    ConstVectorRef angle_axis_ref(angle_axis, 3);
    if (use_quaternions_)
    {
        AngleAxisToQuaternion(angle_axis, camera);
    }
    else
    {
        VectorRef(camera, 3) = angle_axis_ref;
    }

    // t = -R * c
    AngleAxisRotatePoint(angle_axis, center, camera + camera_block_size() - 6);
    VectorRef(camera + camera_block_size() - 6, 3) *= -1.0;
}

void BALProblem::Normalize()
{
    // Compute the marginal median of the geometry
    // 计算几何边缘的中位数
    std::vector<double> tmp(num_points_);
    // 保存点云的中点坐标
    Eigen::Vector3d median;

    // 指向空间三维点
    double *points = mutable_points();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < num_points_; ++j)
        {
            tmp[j] = points[3 * j + i];
        }

        // 分别计算x,y,z的中值
        median(i) = Median(&tmp);
    }

    for (int i = 0; i < num_points_; ++i)
    {
        VectorRef point(points + 3 * i, 3);
        Eigen::Vector3d p = point;
        // 1 范数，代表所有元素的绝对值之和
        tmp[i] = (point - median).lpNorm<1>();
    }

    // 中值绝对值偏差
    const double median_absolute_deviation = Median(&tmp);

    // Scale so that the median absolute deviation of the resulting
    // reconstruction is 100

    const double scale = 100.0 / median_absolute_deviation;

    // X = scale * (X - median)
    for (int i = 0; i < num_points_; ++i)
    {
        VectorRef point(points + 3 * i, 3);
        point = scale * (point - median);
    }

    double *cameras = mutable_cameras();
    double angle_axis[3];
    double center[3];
    for (int i = 0; i < num_cameras_; ++i)
    {
        double *camera = cameras + camera_block_size() * i;
        CameraToAngelAxisAndCenter(camera, angle_axis, center);
        // center = scale * (center - median)
        VectorRef(center, 3) = scale * (VectorRef(center, 3) - median);
        AngleAxisAndCenterToCamera(angle_axis, center, camera);
    }
}

void BALProblem::Perturb(const double rotation_sigma,
                         const double translation_sigma,
                         const double point_sigma)
{
    assert(point_sigma >= 0.0);
    assert(rotation_sigma >= 0.0);
    assert(translation_sigma >= 0.0);

    // 指向空间的三维点
    double *points = mutable_points();
    if (point_sigma > 0)
    {
        for (int i = 0; i < num_points_; ++i)
        {
            PerturbPoint3(point_sigma, points + 3 * i);
        }
    }

    for (int i = 0; i < num_cameras_; ++i)
    {
        double *camera = mutable_cameras() + camera_block_size() * i;

        double angle_axis[3];
        double center[3];
        // Perturb in the rotation of the camera in the angle-axis
        // representation
        CameraToAngelAxisAndCenter(camera, angle_axis, center);
        if (rotation_sigma > 0.0)
        {
            PerturbPoint3(rotation_sigma, angle_axis);
        }
        AngleAxisAndCenterToCamera(angle_axis, center, camera);

        if (translation_sigma > 0.0)
            PerturbPoint3(translation_sigma, camera + camera_block_size() - 6);
    }
}