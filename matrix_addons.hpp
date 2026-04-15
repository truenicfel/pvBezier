
// Synonym for computing the norm
inline RealScalar length() const { return norm(); }
// Synonym for computing the squared norm
inline RealScalar lengthSquared() const { return squaredNorm(); }

// Number of elements
static const int scalar_count = SizeAtCompileTime;
typedef Scalar scalar_type;

// Gets the xy components of the vector
inline Matrix<Scalar, 2, 1> xy() const { return Matrix<Scalar, 2, 1>(this->x(), this->y()); }
// Gets the xyz components of the vector
inline Matrix<Scalar, 3, 1> xyz() const { return Matrix<Scalar, 3, 1>(this->x(), this->y(), this->z()); }

// Returns a pointer to the underlying data
inline const Scalar* ptr() const { return derived().data(); }
// Returns a pointer to the underlying data
inline Scalar* ptr() { return derived().data(); }

// Sorts all entries in ascending order
void sortAscend()
{
    std::sort(derived().data(), derived().data() + size());
}
//
// Sorts all entries in descending order
void sortDescend()
{
    std::sort(derived().data(), derived().data() + size(), std::greater<Scalar>());
}

// Component-wise division
template <typename OtherDerived>
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_quotient_op<Scalar>, const Derived, const OtherDerived>
operator/(const MatrixBase<OtherDerived>& other) const
{
    return cwiseQuotient(other);
}

// Component-wise division
template <typename OtherDerived>
EIGEN_STRONG_INLINE bool
operator<(const EIGEN_CURRENT_STORAGE_BASE_CLASS<OtherDerived>& other) const
{
    for (int i = 0; i < scalar_count; ++i)
        if (this->ptr()[i] != other.ptr()[i])
            return this->ptr()[i] < other.ptr()[i];
    return false;
}

// Computes the dot product.
template <typename OtherDerived>
static Scalar dot(const MatrixBase<OtherDerived>& otherA, const MatrixBase<OtherDerived>& otherB)
{
    return otherA.dot(otherB);
}

// Computes the cross product.
template <typename OtherDerived>
static Matrix<Scalar, RowsAtCompileTime, 1> cross(const MatrixBase<OtherDerived>& otherA, const MatrixBase<OtherDerived>& otherB)
{
    return otherA.cross(otherB);
}

// Linear interpolation
template <typename OtherDerived>
static Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> lerp(const MatrixBase<OtherDerived>& otherA, const MatrixBase<OtherDerived>& otherB, Scalar t)
{
    return Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>(otherA * (1 - t) + otherB * t);
}

// Computes an exponential map
// Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> exponentialMap(const RealScalar dt) const
//{
//	typedef Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> Mat;
//	typedef Matrix<std::complex<Scalar>, RowsAtCompileTime, ColsAtCompileTime> MatC;
//	// compute eigenvalues and eigenvectors
//	Eigen::EigenSolver<Mat> solver(this->derived());
//
//	// place eigenvectors as columns in a matrix
//	MatC V;
//	for (int c = 0; c < ColsAtCompileTime; ++c)
//		for (int r = 0; r < RowsAtCompileTime; ++r)
//			V(r, c) = solver.eigenvectors().col(c)[r];
//
//	// compute the inverse of the eigenvector matrix
//	MatC Vinv = FullPivLU<MatC>(V).inverse();
//
//	// compute matrix exponential on the eigenvalues
//	MatC L;
//	for (int c = 0; c < ColsAtCompileTime; ++c)
//		L(c, c) = std::exp(solver.eigenvalues()[c] * dt);
//
//	// multiply, get real parts and convert to our format
//	MatC N = V * L * Vinv;
//	Mat Nreal;
//	for (int c = 0; c < ColsAtCompileTime; ++c)
//		for (int r = 0; r < RowsAtCompileTime; ++r)
//			Nreal(r, c) = N(r, c).real();
//	return Nreal;
//}

// Linear interpolation
static Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> lerp(const MatrixBase<Derived>& otherA, const MatrixBase<Derived>& otherB, const MatrixBase<Derived>& t)
{
    Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> result;
    for (int c = 0; c < ColsAtCompileTime; ++c)
        for (int r = 0; r < RowsAtCompileTime; ++r)
            result(r, c) = otherA(r, c) * (Scalar(1) - t(r, c)) + otherB(r, c) * t(r, c);
    return result;
}

// Computes the symmetric part of the matrix S = (M + M.transpose()) / 2.
Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> symmetricPart()
{
    return Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>((this->derived() + this->transpose()) / 2);
}

// Computes the antisymmetric part of the matrix S = (M - M.transpose()) / 2.
Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> antiSymmetricPart()
{
    return Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>((this->derived() - this->transpose()) / 2);
}

// Gets a left-handed view matrix from position "eye", in direction "at" with a given "up" vector.
static Matrix<Scalar, 4, 4> lookAtLH(const Matrix<Scalar, 3, 1>& eye, const Matrix<Scalar, 3, 1>& at, const Matrix<Scalar, 3, 1>& up)
{
    Matrix<Scalar, 3, 1> zaxis = (at - eye).normalized();
    Matrix<Scalar, 3, 1> xaxis = up.cross(zaxis).normalized();
    Matrix<Scalar, 3, 1> yaxis = zaxis.cross(xaxis);

    // look at view
    Scalar la_view[] = {
        xaxis.x(), yaxis.x(), zaxis.x(), 0,
        xaxis.y(), yaxis.y(), zaxis.y(), 0,
        xaxis.z(), yaxis.z(), zaxis.z(), 0,
        -xaxis.dot(eye), -yaxis.dot(eye), -zaxis.dot(eye), 1
    };
    return Matrix<Scalar, 4, 4>(la_view);
}

// Gets a right-handed view matrix from position "eye", in direction "at" with a given "up" vector.
static Matrix<Scalar, 4, 4> lookAtRH(const Matrix<Scalar, 3, 1>& eye, const Matrix<Scalar, 3, 1>& at, const Matrix<Scalar, 3, 1>& up)
{
    Matrix<Scalar, 3, 1> zaxis = (eye - at).normalized();
    Matrix<Scalar, 3, 1> xaxis = Matrix<Scalar, 3, 1>::cross(up, zaxis).normalized();
    Matrix<Scalar, 3, 1> yaxis = Matrix<Scalar, 3, 1>::cross(zaxis, xaxis);

    // look at view
    Scalar la_view[] = {
        xaxis.x(), yaxis.x(), zaxis.x(), 0,
        xaxis.y(), yaxis.y(), zaxis.y(), 0,
        xaxis.z(), yaxis.z(), zaxis.z(), 0,
        xaxis.dot(eye), yaxis.dot(eye), zaxis.dot(eye), 1
    };
    return Matrix<Scalar, 4, 4>(la_view);
}

// Gets a left-handed projection matrix from field of view "fovy", aspect ratio, near plane distance "znear" and far plane distance "zfar"
static Matrix<Scalar, 4, 4> perspectiveFovLH(Scalar fovy, Scalar aspect_ratio, Scalar znear, Scalar zfar)
{
    Scalar yScale = 1 / std::tan(fovy / 2);
    Scalar xScale = yScale / aspect_ratio;

    // look at view
    Scalar la_view[] = {
        xScale, 0, 0, 0,
        0, yScale, 0, 0,
        0, 0, zfar / (zfar - znear), 1,
        0, 0, -znear * zfar / (zfar - znear), 0
    };
    return Matrix<Scalar, 4, 4>(la_view);
}

// Gets a right-handed projection matrix from field of view "fovy", aspect ratio, near plane distance "znear" and far plane distance "zfar"
static Matrix<Scalar, 4, 4> perspectiveFovRH(Scalar fovy, Scalar aspect_ratio, Scalar znear, Scalar zfar)
{
    Scalar yScale = 1 / tan(fovy / 2);
    Scalar xScale = yScale / aspect_ratio;

    // look at view
    Scalar la_view[] = {
        xScale, 0, 0, 0,
        0, yScale, 0, 0,
        0, 0, zfar / (znear - zfar), -1,
        0, 0, znear * zfar / (znear - zfar), 0
    };
    return Matrix<Scalar, 4, 4>(la_view);
}

// Gets a left-handed, orthographic projection matrix.
static Matrix<Scalar, 4, 4> orthoLH(Scalar width, Scalar height, Scalar znear, Scalar zfar)
{
    Scalar halfWidth  = width / 2;
    Scalar halfHeight = height / 2;
    return orthoOffCenterLH(-halfWidth, halfWidth, -halfHeight, halfHeight, znear, zfar);
}

// Gets a right-handed, orthographic projection matrix.
static Matrix<Scalar, 4, 4> orthoRH(Scalar width, Scalar height, Scalar znear, Scalar zfar)
{
    Scalar halfWidth  = width / 2;
    Scalar halfHeight = height / 2;
    return orthoOffCenterRH(-halfWidth, halfWidth, -halfHeight, halfHeight, znear, zfar);
}

// Gets a left-handed, customized orthographic projection matrix.
static Matrix<Scalar, 4, 4> orthoOffCenterLH(Scalar left, Scalar right, Scalar bottom, Scalar top, Scalar znear, Scalar zfar)
{
    Scalar zRange    = 1 / (zfar - znear);
    Scalar la_view[] = {
        2 / (right - left), 0, 0, (left + right) / (left - right),
        0, 2 / (top - bottom), 0, (top + bottom) / (bottom - top),
        0, 0, zRange, -znear * zRange,
        0, 0, 0, 1
    };
    return Matrix<Scalar, 4, 4>(la_view);
}

// Gets a right-handed, customized orthographic projection matrix.
static Matrix<Scalar, 4, 4> orthoOffCenterRH(Scalar left, Scalar right, Scalar bottom, Scalar top, Scalar znear, Scalar zfar)
{
    Scalar zRange    = 1 / (zfar - znear);
    Scalar la_view[] = {
        2 / (right - left), 0, 0, (left + right) / (left - right),
        0, 2 / (top - bottom), 0, (top + bottom) / (bottom - top),
        0, 0, -zRange, -znear * zRange,
        0, 0, 0, 1
    };
    return Matrix<Scalar, 4, 4>(la_view);
}

// Gets a translation matrix.
static Matrix<Scalar, 4, 4> translate(const Matrix<Scalar, 3, 1>& translation)
{
    Scalar matrix[] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        translation.x(), translation.y(), translation.z(), 1
    };
    return Matrix<Scalar, 4, 4>(matrix);
}

// Gets a scaling matrix.
static Matrix<Scalar, 4, 4> scale(const Matrix<Scalar, 3, 1>& scaling)
{
    Scalar matrix[] = {
        scaling.x(), 0, 0, 0,
        0, scaling.y(), 0, 0,
        0, 0, scaling.z(), 0,
        0, 0, 0, 1
    };
    return Matrix<Scalar, 4, 4>(matrix);
}

// Gets a diagonal matrix.
static Matrix<Scalar, 4, 4> diag(const Matrix<Scalar, 3, 1>& diagonal)
{
    Scalar matrix[] = {
        diagonal.x(), 0, 0, 0,
        0, diagonal.y(), 0, 0,
        0, 0, diagonal.z(), 0,
        0, 0, 0, 0
    };
    return Matrix<Scalar, 4, 4>(matrix);
}

// Create a perspective transformation. (Maps [near, far] to [0, 1])
// Projects vectors in camera space onto a plane at z=1:
//   x_proj = x / z
//   y_proj = y / z
//   z_proj = (far * (z - near)) / (z * (far-near))
static Matrix<Scalar, 4, 4> perspective(Scalar fov, Scalar near_, Scalar far_)
{
    Scalar recip    = 1.f / (far_ - near_);
    Scalar cot      = Scalar(1) / std::tan((fov * .5f) / 180. * EIGEN_PI);
    Scalar matrix[] = {
        cot, 0, 0, 0,
        0, cot, 0, 0,
        0, 0, far_ * recip, 1,
        0, 0, -near_ * far_ * recip, 0
    };
    return Matrix<Scalar, 4, 4>(matrix);
}

// Gets the translation component, assuming this is a 4x4 matrix.
Matrix<Scalar, 3, 1> translation() const
{
    assert(RowsAtCompileTime == 4 && ColsAtCompileTime == 4);
    return derived().block(0, 3, 3, 1);
}

// Creates a 2x2 matrix
static Matrix<Scalar, 2, 2> make(
    const Scalar& m00, const Scalar& m01,
    const Scalar& m10, const Scalar& m11)
{
    Matrix<Scalar, 2, 2> m;
    m << m00, m01, m10, m11;
    return m;
}

// Creates a 2x2 matrix from column vectors.
static Matrix<Scalar, 2, 2> make(
    const Matrix<Scalar, 2, 1>& c0, const Matrix<Scalar, 2, 1>& c1)
{
    Matrix<Scalar, 2, 2> m;
    m << c0(0), c1(0), c0(1), c1(1);
    return m;
}

// Creates a 3x3 matrix
static Matrix<Scalar, 3, 3> make(
    const Scalar& m00, const Scalar& m01, const Scalar& m02,
    const Scalar& m10, const Scalar& m11, const Scalar& m12,
    const Scalar& m20, const Scalar& m21, const Scalar& m22)
{
    Matrix<Scalar, 3, 3> m;
    m << m00, m01, m02, m10, m11, m12, m20, m21, m22;
    return m;
}

// Creates a 3x3 matrix from column vectors.
static Matrix<Scalar, 3, 3> make(
    const Matrix<Scalar, 3, 1>& c0, const Matrix<Scalar, 3, 1>& c1, const Matrix<Scalar, 3, 1>& c2)
{
    Matrix<Scalar, 3, 3> m;
    m << c0(0), c1(0), c2(0), c0(1), c1(1), c2(1), c0(2), c1(2), c2(2);
    return m;
}

// Creates a 4x4 matrix
static Matrix<Scalar, 4, 4> make(
    const Scalar& m00, const Scalar& m01, const Scalar& m02, const Scalar& m03,
    const Scalar& m10, const Scalar& m11, const Scalar& m12, const Scalar& m13,
    const Scalar& m20, const Scalar& m21, const Scalar& m22, const Scalar& m23,
    const Scalar& m30, const Scalar& m31, const Scalar& m32, const Scalar& m33)
{
    Matrix<Scalar, 4, 4> m;
    m << m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33;
    return m;
}

// Creates a Vector3 by concatenating a Vector2 and a double
static Matrix<Scalar, 3, 1> make(
    const Matrix<Scalar, 2, 1>& xy, const Scalar& z)
{
    Matrix<Scalar, 3, 1> m;
    m << xy.x(), xy.y(), z;
    return m;
}

// Creates a Vector4 by concatenating a Vector3 and a double
static Matrix<Scalar, 4, 1> make(
    const Matrix<Scalar, 3, 1>& xyz, const Scalar& w)
{
    Matrix<Scalar, 4, 1> m;
    m << xyz.x(), xyz.y(), xyz.z(), w;
    return m;
}

// Creates a 6D vector
static Matrix<Scalar, 6, 1> make(
    const Scalar& x, const Scalar& y, const Scalar& z, const Scalar& u, const Scalar& v, const Scalar& w)
{
    Matrix<Scalar, 6, 1> m;
    m << x, y, z, u, v, w;
    return m;
}