#pragma once
#include<iostream>
#include<cmath>
#include<array>

template<class Derived, class T, uint32_t N>
class VectorBase {
public:
	// CRTP
	Derived& derived() { return static_cast<Derived&>(*this); }
	const Derived& derived() const { return static_cast<const Derived&>(*this); }

	T get(uint32_t i) const;
	void set(uint32_t i, T value);

	Derived operator+(const Derived& other) const;
	Derived operator-(const Derived& other) const;
	Derived operator*(T scalar) const;

	T dot(const Derived& other) const;
	T magnitude() const;
	Derived normalize() const;

};

template<class T, uint32_t N>
class Vector :public VectorBase<Vector<T, N>, T, N> {
public:
	friend VectorBase<Vector, T, N>;

	Vector() = default;

	Vector(std::initializer_list<T> init);

	Vector(const Vector& other) :data(other.data) {}

	Vector(Vector&& other)noexcept :data(std::move(other.data)) {
		other.data.fill(T{ 0 });
	}

	Vector& operator=(const Vector& other);

	Vector& operator=(Vector&& other)noexcept;

	T& operator[](uint32_t i);

	const T& operator[](uint32_t i)const;

	// 迭代器支持
	auto begin() { return data.begin(); }
	auto end() { return data.end(); }
	auto begin() const { return data.begin(); }
	auto end() const { return data.end(); }
	auto cbegin() const { return data.cbegin(); }
	auto cend() const { return data.cend(); }
	constexpr uint32_t size() const { return N; }
private:
	std::array<T, N> data{};
};
template <>
class Vector<float, 3> : public VectorBase<Vector<float, 3>, float, 3> {
public:
	using VectorBase::operator*;

	Vector() = default;

	Vector(std::initializer_list<float> init);

	Vector(const Vector& other) :data(other.data) {}

	Vector(Vector&& other)noexcept :data(std::move(other.data)) { other.data.fill(float{ 0 }); }

	Vector& operator=(const Vector& other);

	Vector& operator=(Vector&& other)noexcept;

	float& operator[](uint32_t i);
		                                           
	const float& operator[](uint32_t i)const;

	Vector cross(const Vector& other) const;

	auto begin() { return data.begin(); }
	auto end() { return data.end(); }
	auto begin() const { return data.begin(); }
	auto end() const { return data.end(); }
	auto cbegin() const { return data.cbegin(); }
	auto cend() const { return data.cend(); }
	constexpr uint32_t size() const { return 3; }
private:
	std::array<float, 3> data{};
};

template<class Derived, class T, uint32_t N>
T VectorBase<Derived, T, N >::get(uint32_t i) const {
	if (i < N && i >= 0) {
		return static_cast<const Derived&>(*this)[i];
	}
	else
	{
		std::cout << "error:Vector index out of range" << std::endl;
		return T(0);
	}
}

template<class Derived, class T, uint32_t N>
void VectorBase<Derived, T, N >::set(uint32_t i, T value) {
	if (i < N && i >= 0) {
		static_cast<Derived&>(*this)[i] = value;
	}
	else
	{
		std::cout << "error:Vector index out of range" << std::endl;
	}
}

template<class Derived, class T, uint32_t N>
Derived VectorBase<Derived, T, N >::operator+(const Derived& other) const {
	Derived result{};
	for (uint32_t i = 0; i < N; ++i) {
		result[i] = derived()[i] + other[i];
	}
	return result;
}

template<class Derived, class T, uint32_t N>
Derived VectorBase<Derived, T, N >::operator-(const Derived& other) const {
	Derived result{};
	for (uint32_t i = 0; i < N; ++i) {
		result[i] = derived()[i] - other[i];
	}
	return result;
}

template<class Derived, class T, uint32_t N>
Derived VectorBase<Derived, T, N >::operator*(T scalar) const {
	Derived result{};
	for (uint32_t i = 0; i < N; ++i) {
		result[i] = derived()[i] * scalar;
	}
	return result;
}

template<class Derived, class T, uint32_t N>
T VectorBase<Derived, T, N >::dot(const Derived& other) const {
	T sum = 0;
	for (uint32_t i = 0; i < N; ++i) {
		sum += derived()[i] * other[i];
	}
	return sum;
}

template<class Derived, class T, uint32_t N>
T VectorBase<Derived, T, N >::magnitude() const {
	return std::sqrt(dot(derived()));
}

template<class Derived, class T, uint32_t N>
Derived VectorBase<Derived, T, N >::normalize() const {
	T mag = magnitude();
	if (mag == 0) throw std::runtime_error("Cannot normalize zero vector");
	return derived() * (1 / mag);
}

template<class T, uint32_t N>
Vector<T, N>::Vector(std::initializer_list<T> init) {
	//std::cout<<"初始化列表" << std::endl;
	if (init.size() != N)
		throw std::invalid_argument("Initializer size mismatch");
	std::copy(init.begin(), init.end(), data.begin());
}

template<class T, uint32_t N>
Vector<T, N>& Vector<T, N>::operator=(const Vector& other) {
	//std::cout << "赋值运算符" << std::endl;
	if (this != &other) {
		data = other.data;
	}

	return *this;
}

template<class T, uint32_t N>
Vector<T, N>& Vector<T, N>::operator=(Vector&& other)noexcept {
	//std::cout << "移动赋值" << std::endl;
	if (this != &other) {
		data = std::move(other.data);
		other.data.fill(T{ 0 });
	}

	return *this;
}

template<class T, uint32_t N>
T& Vector<T, N>::operator[](uint32_t i)
{
	if (i >= N) throw std::out_of_range("Vector index out of range");
	return data[i];
}

template<class T, uint32_t N>
const T& Vector<T, N>::operator[](uint32_t i)const
{
	if (i >= N) throw std::out_of_range("Vector index out of range");
	return data[i];
}

inline Vector<float, 3>::Vector(std::initializer_list<float> init) {
	//std::cout<<"初始化列表" << std::endl;
	std::copy(init.begin(), init.end(), data.begin());
}

inline Vector<float, 3>& Vector<float, 3>::operator=(const Vector& other) {
	//std::cout << "拷贝赋值" << std::endl;
	if (this != &other) {
		data = other.data;
	}
	return *this;
}

inline Vector<float, 3>& Vector<float, 3>::operator=(Vector&& other)noexcept {
	//std::cout << "移动赋值" << std::endl;
	if (this != &other) {
		data = std::move(other.data);
		other.data.fill(float{ 0 });
	}
	return *this;
}

inline float& Vector<float, 3>::operator[](uint32_t i)
{
	if (i >= 3) throw std::out_of_range("Vector index out of range");
	return data[i];
}

inline const float& Vector<float, 3>::operator[](uint32_t i)const
{
	if (i >= 3) throw std::out_of_range("Vector index out of range");
	return data[i];
}

inline Vector<float, 3> Vector<float, 3>::cross(const Vector& other) const {
	return {
		data[1] * other[2] - data[2] * other[1],
		data[2] * other[0] - data[0] * other[2],
		data[0] * other[1] - data[1] * other[0]
	};
}

//奇异递归模板模式(CRTP)
template<class Derived, class T, uint32_t Rows, uint32_t Cols>
class MatrixBase {
public:
	Derived& derived() { return static_cast<Derived&>(*this); }
	const Derived& derived() const { return static_cast<const Derived&>(*this); }

	//=== 基础运算符 ===//
	Derived operator+(const Derived& other) const {
		Derived result;
		for (uint32_t i = 0; i < Rows; ++i) {
			result[i] = derived()[i] + other[i];
		}
		return result;
	}

	// 矩阵减法
	Derived operator-(const Derived& other) const {
		Derived result;
		for (uint32_t i = 0; i < Rows; ++i) {
			result[i] = derived()[i] - other[i];
		}
		return result;
	}

	// 标量乘法
	Derived operator*(T scalar) const {
		Derived result;
		for (uint32_t i = 0; i < Rows; ++i) {
			result[i] = derived()[i] * scalar;
		}
		return result;
	}

	// 转置矩阵
	Derived transpose() const {
		Derived transposed;
		for (uint32_t i = 0; i < Rows; ++i) {
			for (uint32_t j = 0; j < Cols; ++j) {
				transposed[j][i] = derived()[i][j];
			}
		}
		return transposed;
	}

	// 归一化
	Derived normalize() const {
		Derived result;
		for (uint32_t i = 0; i < Rows; ++i) {
			T length = 0;
			for (uint32_t j = 0; j < Cols; ++j) {
				length += derived()[i][j] * derived()[i][j];
			}
			length = std::sqrt(length);

			if (length == 0) throw std::runtime_error("Zero vector cannot be normalized");

			for (uint32_t j = 0; j < Cols; ++j) {
				result[i][j] = derived()[i][j] / length;
			}
		}
		return result;
	}

	//=== 静态方法 ===//
	static Derived identity()
	{
		static_assert(Rows == Cols, "Identity matrix must be square");
		Derived mat;
		for (uint32_t i = 0; i < Rows; ++i) {
			mat[i][i] = T(1);
		}
		return mat;
	}

	static Derived zero() 
	{
		Derived mat;
		for (uint32_t i = 0; i < Rows; ++i) {
			for (uint32_t j = 0; j < Cols; ++j) {
				mat[i][j] = T(0);
			}
		}
		return mat;
	}
};


/**
	* @brief 矩阵类模板
	* @tparam T    元素类型
	* @tparam Rows 行数
	* @tparam Cols 列数
**/
template <class T, uint32_t Rows, uint32_t Cols>
class Matrix : public MatrixBase<Matrix<T, Rows, Cols>, T, Rows, Cols> {
public:
	//=== 构造函数 ===//
	Matrix() {
		// 初始化为零矩阵
		for (uint32_t row = 0; row < Rows; ++row) {
			mMatrix[row] = Vector<T, Cols>(); // 每行有Cols列
		}
		// 方阵则初始化为单位矩阵
		if constexpr (Rows == Cols) {
			for (uint32_t i = 0; i < Rows; ++i) {
				mMatrix[i][i] = T(1);
			}
		}
	}

	// 可变参数模板构造函数（排除拷贝场景）
	template <
		typename... Vectors,
		typename = std::enable_if_t<(!std::is_same_v<std::decay_t<Vectors>, Matrix> && ...)>>
		explicit Matrix(Vectors&&... vectors) {
		static_assert(sizeof...(Vectors) == Rows,
			"Number of vectors must match matrix rows (Rows)");
		initMatrix(0, try_create_vector(std::forward<Vectors>(vectors))...);
	}

	// 初始化列表构造
	Matrix(std::initializer_list<std::initializer_list<T>> init) {
		if (init.size() != Rows) {
			throw std::invalid_argument("Row count mismatch");
		}
		auto rowIt = init.begin();
		for (uint32_t i = 0; i < Rows; ++i, ++rowIt) {
			if (rowIt->size() != Cols) {
				throw std::invalid_argument("Column count mismatch in row " + std::to_string(i));
			}
			std::copy(rowIt->begin(), rowIt->end(), mMatrix[i].begin());
		}
	}

	//拷贝构造
	Matrix(const Matrix<T, Rows, Cols>& matrix) {
		for (int i = 0; i < Rows; ++i)
		{
			mMatrix[i] = matrix[i];
		}
	}

	//移动构造
	Matrix(Matrix<T, Rows, Cols>&& matrix) noexcept {
		for (uint32_t i = 0; i < Rows; ++i) {
			mMatrix[i] = std::move(matrix[i]);
		}
	}

	//=== 基本操作 ===//
	void setRow(uint32_t row, const Vector<T, Cols>& vector) {
		if (row >= Rows) { // 检查行索引范围 [0, Rows)
			throw std::out_of_range("Row index out of bounds");
		}
		mMatrix[row] = vector;
	}

	void setColumn(uint32_t column, const Vector<T, Rows>& vector) {
		if (column >= Cols) { // 检查列索引范围 [0, Cols)
			throw std::out_of_range("Column index out of bounds");
		}
		for (uint32_t i = 0; i < Rows; ++i) {
			mMatrix[i][column] = vector[i];
		}
	}

	Vector<T, Cols>& getRow(uint32_t row) {
		if (row >= Rows) {
			throw std::out_of_range("Row index out of bounds");
		}
		return mMatrix[row];
	}

	Vector<T, Rows> getColumn(uint32_t column) const {
		if (column >= Cols) {
			throw std::out_of_range("Column index out of bounds");
		}
		Vector<T, Rows> col;
		for (uint32_t i = 0; i < Rows; ++i) {
			col[i] = mMatrix[i][column];
		}
		return col;
	}

	//=== 运算符重载 ===//
	Matrix& operator=(const Matrix& other) noexcept {
		if (this != &other) {
			for (uint32_t i = 0; i < Rows; ++i) { // 循环Rows次
				mMatrix[i] = other[i];
			}
		}
		return *this;
	}

	Matrix& operator=(Matrix&& other) noexcept {
		for (uint32_t i = 0; i < Rows; ++i) { // 循环Rows次
			mMatrix[i] = std::move(other[i]);
		}
		return *this;
	}

	// 行访问（返回Cols列的向量）
	Vector<T, Cols>& operator[](uint32_t row) {
		return mMatrix[row];
	}
	const Vector<T, Cols>& operator[](uint32_t row) const {
		return mMatrix[row];
	}

	auto begin() { return mMatrix.begin(); }
	auto end() { return mMatrix.end(); }
	auto begin() const { return mMatrix.begin(); }
	auto end() const { return mMatrix.end(); }

private:
	//=== 辅助函数 ===//
	template <typename U>
	Vector<T, Cols> try_create_vector(U&& arg) {
		static_assert(
			std::is_constructible_v<Vector<T, Cols>, U>,
			"Argument must be convertible to Vector<T, Cols>"
			);
		return Vector<T, Cols>(std::forward<U>(arg));
	}

	template <typename First, typename... Rest>
	void initMatrix(uint32_t index, First&& first, Rest&&... rest) {
		mMatrix[index] = std::forward<First>(first);
		if constexpr (sizeof...(Rest) > 0) {
			initMatrix(index + 1, std::forward<Rest>(rest)...);
		}
	}

	void initMatrix(uint32_t) {} // 终止递归

private:
	// 存储结构：Rows行，每行是Cols列的向量
	Vector<Vector<T, Cols>, Rows> mMatrix;
};

// 矩阵乘法（通用实现）
template <typename T, uint32_t Rows, uint32_t N, uint32_t Cols>
Matrix<T, Rows, Cols> operator*(
	const Matrix<T, Rows, N>& lhs,
	const Matrix<T, N, Cols>& rhs) {
	Matrix<T, Rows, Cols> result;
	for (uint32_t i = 0; i < Rows; ++i) {
		for (uint32_t j = 0; j < Cols; ++j) {
			T sum = 0;
			for (uint32_t k = 0; k < N; ++k) {
				sum += lhs[i][k] * rhs[k][j];
			}
			result[i][j] = sum;
		}
	}
	return result;
}

//逆矩阵（仅适用于方阵）
template <class T, uint32_t N>
Matrix<T, N, N> inverse(const Matrix<T, N, N>& mat) {
	const T epsilon = std::numeric_limits<T>::epsilon();

	// 创建增广矩阵 [A|I]
	Matrix<T, N, 2 * N> aug;

	// 填充原矩阵部分
	for (uint32_t i = 0; i < N; ++i) {
		for (uint32_t j = 0; j < N; ++j) {
			aug[i][j] = mat[i][j];
		}
		aug[i][i + N] = T(1); // 单位矩阵部分
	}

	// 高斯-约旦消元
	for (uint32_t pivot = 0; pivot < N; ++pivot) {
		// 寻找主元行
		uint32_t max_row = pivot;
		T max_val = std::abs(aug[pivot][pivot]);
		for (uint32_t row = pivot + 1; row < N; ++row) {
			T current = std::abs(aug[row][pivot]);
			if (current > max_val) {
				max_val = current;
				max_row = row;
			}
		}

		// 奇异矩阵检测
		if (max_val < epsilon) {
			throw std::runtime_error("Matrix is singular");
		}

		// 行交换
		if (max_row != pivot) {
			using std::swap;
			swap(aug[pivot], aug[max_row]);
		}

		// 归一化主元行
		T diag = aug[pivot][pivot];
		for (uint32_t col = pivot; col < 2 * N; ++col) {
			aug[pivot][col] /= diag;
		}

		// 消元其他行
		for (uint32_t row = 0; row < N; ++row) {
			if (row != pivot && std::abs(aug[row][pivot]) > epsilon) {
				T factor = aug[row][pivot];
				for (uint32_t col = pivot; col < 2 * N; ++col) {
					aug[row][col] -= factor * aug[pivot][col];
				}
			}
		}
	}

	// 提取逆矩阵部分
	Matrix<T, N, N> inv;
	for (uint32_t i = 0; i < N; ++i) {
		for (uint32_t j = 0; j < N; ++j) {
			inv[i][j] = aug[i][j + N];
		}
	}

	return inv;
}

template<class T, uint32_t Rows, uint32_t Cols>
void print(Matrix<T, Rows, Cols>& matrix) {
	for (uint32_t i = 0; i < Rows; i++)
	{
		for (uint32_t j = 0; j < Cols; j++)
		{
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

using float2 = Vector<float, 2>;
using float3 = Vector<float, 3>;
using float4 = Vector<float, 4>;
using int2 = Vector<int, 2>;
using int3 = Vector<int, 3>;
using int4 = Vector<int, 4>;
using Matrixf4x4 = Matrix<float, 4, 4>;
using Matrixf3x3 = Matrix<float, 3, 3>;

inline int cross(int2& v1, int2& v2) {

	return v1[0] * v2[1] - v1[1] * v2[0];
}

inline float cross(float2& v1, float2& v2) {

	return v1[0] * v2[1] - v1[1] * v2[0];
}

struct screenPoint {
	int2 screenVertex{ 0,0 };
	int4 screenColor{ 0,0,0,0 };
};

class Vertex {
public:
	Vertex() {
		for (uint32_t i = 0; i < 2; i++)
		{
			vertexCoord[i] = 0;
		}
		for (uint32_t i = 0; i < 3; i++)
		{
			vertexPosition[i] = 0;
		}
		for (uint32_t i = 0; i < 4; i++)
		{
			vertexColor[i] = 0;
		}
	}
	Vertex(Vertex& vertex) {
		vertexPosition = vertex.getPosition();
		vertexCoord = vertex.getVertexCoord();
		vertexColor = vertex.getColor();
	}
	Vertex(const Vertex& vertex) {
		vertexPosition = vertex.getPosition();
		vertexCoord = vertex.getVertexCoord();
		vertexColor = vertex.getColor();
	}

	Vertex& setColor(float x, float y, float z, float w) {
		vertexColor[0] = x;
		vertexColor[1] = y;
		vertexColor[2] = z;
		vertexColor[3] = w;
		return *this; // 返回当前对象的引用
	}

	Vertex& setColor(const float4& Color){
		vertexColor[0] = Color[0];
		vertexColor[1] = Color[1];
		vertexColor[2] = Color[2];
		vertexColor[3] = Color[3];
		return *this; // 返回当前对象的引用
	}

	Vertex& setPosition(float x, float y, float z) {
		vertexPosition[0] = x;
		vertexPosition[1] = y;
		vertexPosition[2] = z;
		return *this;
	}

	Vertex& setVertexCoord(float x, float y) {
		vertexCoord[0] = x;
		vertexCoord[1] = y;
		return *this;
	}
	float4 getColor()const {
		return vertexColor;
	}
	float3 getPosition()const {
		return vertexPosition;
	}
	float2 getVertexCoord()const {
		return vertexCoord;
	}
	void setPosition(float3& position) {
		vertexPosition = position;
	}
	void setVertexCoord(float2& VertexCoord) {
		vertexCoord = VertexCoord;
	}

	~Vertex() {}
private:
	float2  vertexCoord;
	float3  vertexPosition;
	float4  vertexColor;
};

template <class T>
Vector<T, 4> operator*(const Matrix<T, 4, 4>& mat, const Vector<T, 3>& vec) {
	// 转换为齐次坐标（w=1）
	Vector<T, 4> pos = Vector<T, 4>{ vec[0],vec[1],vec[2],1 };

	Vector<T, 4> result;
	for (uint32_t row = 0; row < 4; ++row) {
		result[row] = mat[row].dot( pos);
	}
	return result;
}

// 四元数
template<typename T>
class Quaternion {
public:
	T x, y, z, w;

	Matrix<T, 4, 4> ToRotationMatrix() const {
		Matrix<T, 4, 4> mat;
		T xx = x * x, yy = y * y, zz = z * z;
		T xy = x * y, xz = x * z, yz = y * z;
		T wx = w * x, wy = w * y, wz = w * z;

		mat[0][0] = 1 - 2 * (yy + zz);
		mat[0][1] = 2 * (xy - wz);
		mat[0][2] = 2 * (xz + wy);

		mat[1][0] = 2 * (xy + wz);
		mat[1][1] = 1 - 2 * (xx + zz);
		mat[1][2] = 2 * (yz - wx);

		mat[2][0] = 2 * (xz - wy);
		mat[2][1] = 2 * (yz + wx);
		mat[2][2] = 1 - 2 * (xx + yy);

		mat[3][3] = 1;
		return mat;
	}
};

// 欧拉角转四元数
template<typename T>
Quaternion<T> EulerToQuaternion(const Vector<T, 3>& euler) {
	T cy = cos(euler[2] * 0.5); // Yaw
	T sy = sin(euler[2] * 0.5);
	T cp = cos(euler[1] * 0.5); // Pitch
	T sp = sin(euler[1] * 0.5);
	T cr = cos(euler[0] * 0.5); // Roll
	T sr = sin(euler[0] * 0.5);

	Quaternion<T> q;
	q.w = cr * cp * cy + sr * sp * sy;
	q.x = sr * cp * cy - cr * sp * sy;
	q.y = cr * sp * cy + sr * cp * sy;
	q.z = cr * cp * sy - sr * sp * cy;
	return q;
}

template<typename T>
Matrix<T, 4, 4> CreateModelMatrix(
	const Vector<T, 3>& translation,
	const Vector<T, 3>& rotationEuler, // 欧拉角（弧度）XYZ顺序
	const Vector<T, 3>& scale)
{
	// 创建单位矩阵
	Matrix<T, 4, 4> model = Matrix<T, 4, 4>::identity();

	// 应用缩放
	Matrix<T, 4, 4> scaleMat = Matrix<T, 4, 4>::identity();
	scaleMat[0][0] = scale[0];
	scaleMat[1][1] = scale[1];
	scaleMat[2][2] = scale[2];
	model = scaleMat * model;

	// 应用旋转（欧拉角转四元数）
	Quaternion<T> q = EulerToQuaternion(rotationEuler);
	model = q.ToRotationMatrix() * model;

	// 应用平移
	Matrix<T, 4, 4> transMat = Matrix<T, 4, 4>::identity();
	transMat[0][3] = translation[0];
	transMat[1][3] = translation[1];
	transMat[2][3] = translation[2];
	model = transMat * model;

	return model;
}

template<typename T>
Matrix<T, 4, 4> CreateViewMatrix(
	const Vector<T, 3>& eyePos,    // 摄像机位置
	const Vector<T, 3>& target,    // 观察目标点
	const Vector<T, 3>& upWorld)   // 世界空间上方向
{
	Vector<T, 3> forward = (target - eyePos).normalize();
	Vector<T, 3> right = forward.cross(upWorld).normalize();
	Vector<T, 3> up = right.cross(forward);

	Matrix<T, 4, 4> view;
	// 旋转部分
	view[0][0] = right[0];   view[0][1] = right[1];   view[0][2] = right[2];
	view[1][0] = up[0];      view[1][1] = up[1];      view[1][2] = up[2];
	view[2][0] = -forward[0]; view[2][1] = -forward[1]; view[2][2] = -forward[2];

	// 平移部分
	view[0][3] = -right.dot(eyePos);
	view[1][3] = -up.dot(eyePos);
	view[2][3] = forward.dot(eyePos);

	// 齐次坐标
	view[3][3] = 1.0;

	return view;
}

template<typename T>
Matrix<T, 4, 4> CreatePerspectiveProjection(
	T fovY,         // 垂直视野角（弧度）
	T aspectRatio,  // 宽高比（width/height）
	T nearZ,        // 近裁剪面
	T farZ)         // 远裁剪面
{
	Matrix<T, 4, 4> proj;
	T tanHalfFov = tan(fovY / 2);

	proj[0][0] = 1 / (aspectRatio * tanHalfFov);
	proj[1][1] = 1 / tanHalfFov;
	proj[2][2] = -(farZ + nearZ) / (farZ - nearZ);
	proj[2][3] = -2 * farZ * nearZ / (farZ - nearZ);
	proj[3][2] = -1; 
	proj[3][3] = 0;

	return proj;
}

template<typename T>
Matrix<T, 4, 4> CreateOrthographicProjection(
	T left, T right,
	T bottom, T top,
	T nearZ, T farZ)
{
	Matrix<T, 4, 4> proj;

	proj[0][0] = 2 / (right - left);
	proj[1][1] = 2 / (top - bottom);
	proj[2][2] = -2 / (farZ - nearZ);

	proj[0][3] = -(right + left) / (right - left);
	proj[1][3] = -(top + bottom) / (top - bottom);
	proj[2][3] = -(farZ + nearZ) / (farZ - nearZ);
	proj[3][3] = 1;

	return proj;
}

inline float edgeFunction(const float2& a, const float2& b, const float2& p) {
	return (p[0] - a[0]) * (b[1] - a[1]) - (p[1] - a[1]) * (b[0] - a[0]);
}

struct VertexShaderOutPut
{
	float4	clipPos;
	float4  vertexColor;
	float2  vertexCoord;
};

struct FragmentShaderOutPut {
	int4 pixelColor;
	int2 screenPosition;
};

struct VertexShaderToRater {
	float4 pixelColor;
	float2 screenPosition;
};



