#ifndef MATRIX_ALG_HPP
#define MATRIX_ALG_HPP

#include <iostream>
#include <vector>
#include "constants.hpp"

struct AppendList{
    int i;
    int j;
    double value;
    AppendList *next;
    AppendList(){i=0; j=0; value = 0.0; next = nullptr;}
};
// /* -------------------------------------------------------------------------------------------------------------------- */
// /* -------- Class for defining an Array of classes  ------------------------------------------------------------------- */
// /* -------------------------------------------------------------------------------------------------------------------- */

// template<typename T, const int N>
// class Array{
// private:
//     T m_Array[N];
// public:
//     Array(); // default constructor
//     void init_();
// };

// template<typename T, const int N>
// Array<T,N>::Array(){
//     init_();
// }

// template<typename T, const int N>
// void Array<T,N>::init_(){
//     for (int i=0; i<N; i++){
//         m_Array[i](); // assuming array of these classes have default constructors
//     }
// }

/* -------------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------------- */
template<typename T, const int N>
class Vector{
    private:
        T m_Vector[N];
    public:
        Vector(); // Constructor
        Vector(T *v, int Size); // Constructor to initialize with values
        void init_();
        int getLength_(){
            return N;
        }
        void setValue(int idx, T value);
        void setVector(T *vec, int Size);
        T getValue(int idx);
        void Add(T val, Vector<T,N> &v_out);
        void Add(Vector<T,N> v_in, Vector<T,N> &v_out);
        void Subtract(T val, Vector<T,N> &v_out);
        void Subtract(Vector<T,N> v_in, Vector<T,N> &v_out);
        void Scale(T scale, Vector<T,N> &v_out);
        void ElementMultiplication(Vector<T,N> v_in, Vector<T,N> &v_out);
        T dotProduct(Vector<T,N> a);
        T *getData();
        void AssignFromList(AppendList *head);
        void displayVector();

};

/// default constructor
template<typename T, const int N>
Vector<T,N>::Vector(){
   init_();
}

/// Initialize the Vector object with a pointer to an array
template<typename T, const int N>
Vector<T,N>::Vector(T *v, int Size){
    if(N > Size){
        init_();
        std::cerr << "All values wont be initialized since size of vector > selected vector size" << std::endl;
    }
    else{
        for(int i=0; i<N; i++){
            m_Vector[i] = v[i];
        }
    }
}

/// Initialize the Vector object to zeros
template<typename T, const int N>
void Vector<T,N>::init_(){
    for (int i=0; i<N; i++){
        setValue(i,0);
    }
}

/// set Value at index i to (value) of current Vector object
template<typename T, const int N>
void Vector<T,N>::setValue(int idx, T value){
    if(idx >= N){ std::cerr << "Index length greater than N" << std::endl;}
    else{ m_Vector[idx] = value; }
}

/// Set the entire vector in current Vector object to the pointer *vec 
template<typename T, const int N>
void Vector<T,N>::setVector(T *vec, int Size){
    if(N > Size){
        init_();
        std::cerr << "All values wont be initialized since size of vector > selected vector size" << std::endl;
    }
    else{
        for(int i=0; i<N; i++){
            m_Vector[i] = vec[i];
        }
    }
}

/// get Value of certain index of current Vector object
template<typename T, const int N> 
T Vector<T,N>::getValue(int idx){
    if(idx>=N){std::cerr << "Index youre trying to access exceeds vector length" << std::endl; return (T)0.0;}
    else{ return m_Vector[idx]; }
}

/// Adds the scalar value (val) to all elements of current Vector object and assigns it
/// to v_out. v_out = current Vector object + (val)
template<typename T, const int N>
void Vector<T,N>::Add(T val, Vector<T,N> &v_out){
    for(int i=0; i<N; i++){
        v_out.setValue(i,m_Vector[i]+val);
    }
}

/// Adds current vector object to v_in and assigns it to v_out
/// v_out = (current Vector object) + v_in
template<typename T, const int N>
void Vector<T,N>::Add(Vector<T,N> v_in, Vector<T,N> &v_out){
    for(int i=0; i<N; i++){
        v_out.setValue(i,m_Vector[i]+v_in.getValue(i));
    }
}

/// Subtract all values of current Vector object with scalar val and assigns it to 
/// v_out. v_out = (current Vector object) - (scalar val)
template<typename T, const int N>
void Vector<T,N>::Subtract(T val, Vector<T,N> &v_out){
    for(int i=0; i<N; i++){
        v_out.setValue(i,m_Vector[i]-val);
    }
}

/// Subtract the current Vector object with the Vector v_in and assigns it to v_out
/// v_out = (current Vector object) - v_in
template<typename T, const int N>
void Vector<T,N>::Subtract(Vector<T,N> v_in, Vector<T,N> &v_out){
    for(int i=0; i<N; i++){
        v_out.setValue(i,m_Vector[i]-v_in.getValue(i));
    }
}

/// Scales the current Vector object by (scale) and assigns it to Vector v_out
/// v_out = (scale)*(current Vector object)
template<typename T, const int N>
void Vector<T,N>::Scale(T scale, Vector<T,N> &v_out){
    for(int i=0; i<N; i++){
        v_out.setValue(i,m_Vector[i] * scale);
    }
}

/// Does element by element multiplication of the values of the current Vector and
/// Vector v_in and then assigns that value to Vector v_out
template<typename T, const int N>
void Vector<T,N>::ElementMultiplication(Vector<T,N> v_in, Vector<T,N> &v_out){
    for(int i=0; i<N; i++){
        v_out.setValue(i, getValue(i)*v_in.getValue(i)); 
    }
}

/// Dot Product of two vectors is performed.
/// scalar value = (current Vector object)(.) a
template<typename T, const int N>
T Vector<T,N>::dotProduct(Vector<T,N> a){
    T sum= (T)0.0;
    for (int i=0; i<N; i++){
        sum += a.getValue(i)*getValue(i);
    }
    return sum;
}

/// gives the pointer to the place where m_Vector is stored. 
template<typename T, const int N>
T *Vector<T,N>::getData(){
    T *temp_vector = NULL;
    temp_vector = m_Vector;
    return temp_vector;
}

template<typename T, const int N>
void Vector<T,N>::AssignFromList(AppendList * head){
    int i;
    for (;head != nullptr;head=head->next){
        T temp = getValue(head->i);
        i = head->i;
        setValue(head->i, temp+(T)head->value);
    }
}

/// Prints vector
template<typename T, const int N>
void Vector<T,N>::displayVector(){
    for(int i=0; i<N; i++){std::cout<<m_Vector[i]<<" ";}
}
/* ---------------------------------------------------------------------------------------------------------------------  */
/* ---------------------------------------------------------------------------------------------------------------------  */
template<typename T, const int m, const int n>
class Matrix{
    private:
        T m_Matrix[m][n];
        bool SQUARE;
    public:
        Matrix(); // Constructor
        Matrix(T *rowMajor, int Size); // imports matrix data from a rowMajor format
        void getSize_(){
            std::cout << "The size is:" << m << " x " << n << std::endl;
        }
        Matrix(T i); // Matrix to generate i*Identity matrix 
        Matrix(Vector<T,m> v1, Vector<T,n> v2);
        void init_();
        void setValue(int i, int j, T val);
        T getValue(int i, int j);
        void setRow(int row_idx, Vector<T,n> r);
        void setColumn(int col_idx, Vector<T,m> c);
        void getRow(int row_idx, Vector<T,n> &r);
        void getColumn(int col_idx, Vector<T,m> &c);
        void Transpose(Matrix<T,n,m> &m_Matrix_T);
        void Add(T val, Matrix<T,m,n> &M_out);
        void Add(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out);
        void Subtract(T val, Matrix<T,m,n> &M_out);
        void Subtract(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out);
        void Multiply(T scale, Matrix<T,m,n> &M_out);
        void Multiply(Vector<T,n> u, Vector<T,m> &u_out);
        template<const int C>
        void Multiply(Matrix<T,n,C> M_in, Matrix<T,m,C> &M_out);
        void ElementMultiplication(Matrix<T,m,n> m1, Matrix<T,m,n> &m2);
        bool checkSquare(){return SQUARE;}
        void Invert(Matrix<T,n,n> &inv);
        T *exportRowMajor();
        void setRowMajor(T *rowMajor, int Size);
        void AssignFromList(AppendList *head);
        void displayMatrix();
};

/// Constructor to initialize the Matrix object with 0 values at all entries
template<typename T, const int m, const int n>
Matrix<T,m,n>::Matrix(){
    init_();
    
}

/// Constructor to generate a Matrix object from a rowMajor array 
template<typename T, const int m, const int n>
Matrix<T,m,n>::Matrix(T *rowMajor, int Size){
    if (m*n != Size){
        std::cerr<<"Matrix size is incompatible" << std::endl;
    }
    else {
        int idx = 0;
        for (int i=0; i<m; i++){
            for (int j=0; j<n; j++){
                m_Matrix[i][j] = rowMajor[idx];
                ++idx;
                }
        }
        SQUARE = (m==n)?true:false;
    }
}

/// Constructor to generate a Matrix object which is i*Identity matrix
template<typename T, const int m, const int n>
Matrix<T,m,n>::Matrix(T i){
    init_();
    if(SQUARE){
        for(int idx=0; idx<n; idx++){
            m_Matrix[idx][idx] = (i)*(T)(1.0);
        }
    }
    else{
        std::cerr<<"Not a square matrix and hence it is initialized with zeros" << std::endl;
    }
}

/// Constructor for object Matrix where current Matrix = v1*v2'
template<typename T, const int m, const int n>
Matrix<T,m,n>::Matrix(Vector<T,m> v1, Vector<T,n> v2){
    Vector<T,m> temp;
    for (int i=0; i<n; i++){
        v1.Scale(v2.getValue(i),temp);
        setColumn(i,temp);
    }
}

/// Initialize all values with values 0.0
template<typename T, const int m, const int n>
void Matrix<T,m,n>::init_(){
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            m_Matrix[i][j] = (T)(0.0);
        }
    }
    SQUARE = (m==n)?true:false;
}

/// Set one single value of the current Matrix object
template<typename T, const int m, const int n>
void Matrix<T,m,n>::setValue(int i, int j, T val){
    if(i>=m || j>=n){
        std::cerr << "Index exceeds matrix dimensions" << std::endl;
    }
    else {
        m_Matrix[i][j] = val;
    }
}

/// Get value from a index i,j
template<typename T, const int m, const int n>
T Matrix<T,m,n>::getValue(int i, int j){
    if(i>=m || j>=n){
        std::cerr<<"Trying to access invalid index values" << std::endl;
        return (T)0.0;
    }
    else{
        return m_Matrix[i][j];
    }
}

/// Set row of current matrix of index row_idx with Vector r
template<typename T, const int m, const int n>
void Matrix<T,m,n>::setRow(int row_idx, Vector<T,n> r){
    if (row_idx >= m){
        std::cerr << "Row index exceeds Matrix dimensions" << std::endl;
    }
    else {
        for(int i = 0; i < n; i++){
            m_Matrix[row_idx][i] = r.getValue(i);
        }
    }
}

/// Set column of current matrix of index col_idx with Vector c
template<typename T, const int m, const int n>
void Matrix<T,m,n>::setColumn(int col_idx, Vector<T,m> c){
    if (col_idx >= n){
        std::cerr << "Column index exceeds Matrix dimensions" << std::endl;
    }
    else {
        for (int i = 0; i < m; i++){
            m_Matrix[i][col_idx] = c.getValue(i);
        }
    }
}

/// Get a row of current Matrix object and assign it to Vector r
template<typename T, const int m, const int n>
void Matrix<T,m,n>::getRow(int row_idx, Vector<T,n> &r){
    if(row_idx >= m){
        std::cerr << "Row index exceeds matrix dimensions" << std::endl;
    }
    else {
        for(int i=0; i<n; i++){
            r.setValue(i,m_Matrix[row_idx][i]);
        }
    }
}

/// Get a column of current Matrix object and assign it to Vector c
template<typename T, const int m, const int n>
void Matrix<T,m,n>::getColumn(int col_idx, Vector<T,m> &c){
    if(col_idx >= n){
        std::cerr << "Column index exceeds matrix dimensions" << std::endl;
    }
    else {
        for(int i=0; i<m; i++){
            c.setValue(i,m_Matrix[i][col_idx]);
        }
    }
}

/// Transposes current Matrix object and assigns it to m_Matrix_T
/// m_Matrix_T = (current Matrix)'
template<typename T, const int m, const int n>
void Matrix<T,m,n>::Transpose(Matrix<T,n,m> &m_Matrix_T){
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            m_Matrix_T.setValue(j,i,m_Matrix[i][j]);
        }
    }
}

/// Add a scalar value (val) and add it to all elements of current Matrix object
/// and assigns it to M_out. M_out = (current Matrix) + val
template<typename T, const int m, const int n>
void Matrix<T,m,n>::Add(T val, Matrix<T,m,n> &M_out){
    for(int i=0; i<m; i++){
        Vector<T,n> r;
        getRow(i,r);
        r.Add(val,r);
        M_out.setRow(i,r);
    }
}

/// Adds current matrix object to Matrix object M_in and assigns to M_out
/// M_out = (current Matrix) + M_in
template<typename T, const int m, const int n>
void Matrix<T,m,n>::Add(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out){
    T *rm1 = exportRowMajor();
    T *rm2 = M_in.exportRowMajor();
    T *vec = (T*)malloc((m*n)*sizeof(T));
    for (int i=0; i < m*n; i++){
        vec[i] = rm1[i] + rm2[i];
    }
    Matrix<T,m,n> test(vec, m*n);
    //M_out.setRowMajor(vec, m*n);
    M_out = test;
}

/// Subtracts a single value from all elements of the current Matrix object and assigns to 
/// M_out. M_out = M_in - (value)
template<typename T, const int m, const int n>
void Matrix<T,m,n>::Subtract(T val, Matrix<T,m,n> &M_out){
    for(int i=0; i<m; i++){
        Vector<T,n> r;
        getRow(i,r);
        r.Subtract(val,r);
        M_out.setRow(i,r);
    }
}

/// Subtracts the parameter M_in from current Matrix object and assigns it to Matrix object
/// M_out. M_out = (current Matrix) - M_in
template<typename T, const int m, const int n>
void Matrix<T,m,n>::Subtract(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out){
    T *rm1 = exportRowMajor();
    T *rm2 = M_in.exportRowMajor();
    T *vec = (T*)malloc((m*n)*sizeof(T));
    for (int i=0; i < m*n; i++){
        vec[i] = rm1[i] - rm2[i];
    }
    Matrix<T,m,n> test(vec);
    //M_out.setRowMajor(vec, m*n);
    M_out = test;
}

/// Scales each element of current Matrix object by (value scale) and assigns it to
/// Matrix object M_out. M_out = (scale)*(currentMatrix)
template<typename T, const int m, const int n>
void Matrix<T,m,n>::Multiply(T scale, Matrix<T,m,n> &M_out){
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            M_out.setValue(i,j,m_Matrix[i][j]*scale);
        }
    }
}

/// This multiplies current Matrix object with a vector u and assigns it to 
/// another vector u_out. u_out = (currentMatrix)*u 
template<typename T, const int m, const int n>
void Matrix<T,m,n>::Multiply(Vector<T,n> u, Vector<T,m> &u_out){
    for(int i=0; i<m; i++){
        Vector<T,n> row ;
        getRow(i,row);
        T val = u.dotProduct(row);
        u_out.setValue(i,val);
    }
}

/// Current Matrix object is matrix multiplied with Matrix object M_in and assigned to
/// Matrix object M_out. M_out = (current)*M_in 
template<typename T, const int m, const int n>
template<const int C>
void Matrix<T,m,n>::Multiply(Matrix<T,n,C> M_in, Matrix<T,m,C> &M_out){
    for(int i=0; i<m; i++){
        Vector<T,n> row;
        getRow(i,row);
        for(int j=0; j<C; j++){
            Vector<T,n> column;
            M_in.getColumn(j,column);
            T val = row.dotProduct(column);
            M_out.setValue(i,j,val);
        }
    }
}

/// Element by Element multiplication is performed on current Matrix object and Matrix
/// object M_in and assigns it to Matrix object M_out
template<typename T, const int m, const int n>
void Matrix<T,m,n>::ElementMultiplication(Matrix<T,m,n> M_in, Matrix<T,m,n> &M_out){
    Vector<T,m> temp1, temp2, temp3;
    for(int j=0; j<n; j++){
        getColumn(j,temp1);
        M_in.getColumn(j,temp2);
        temp1.ElementMultiplication(temp2,temp3);
        M_out.setColumn(j,temp3);
    }
}

/// Inverts the current matrix and assigns it to Matrix object inverse
template<typename T, const int m, const int n>
void Matrix<T,m,n>::Invert(Matrix<T,n,n> &inv){
    if(SQUARE){

    }
    else {
        std::cerr<< "Matrix is not a square matrix and cannot be inverted" << std::endl;
    }
}

/// exports the matrix as a pointer to an array. The rows of the matrix are organized one next to each other.
template<typename T, const int m, const int n>
T *Matrix<T,m,n>::exportRowMajor(){
    T *rowMajor = (T*)malloc((m*n)*sizeof(T));
    int idx = 0;
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            rowMajor[idx] = m_Matrix[i][j];
            ++idx;
        }
    }
    return rowMajor;
}


/// Rearrange a row major vector into matrix format
template<typename T, const int m, const int n>
void Matrix<T,m,n>::setRowMajor(T *rowMajor, int Size){
    if(m*n != Size){
        std::cerr << "Array size is not compatible with the matrix for allocation" << std::endl;
    }
    else {
        int idx = 0;
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                m_Matrix[i][j] = rowMajor[idx];
                ++idx;
            }
        }
    }
}

/// Using LinkedList to append(+add) values to Matrix elements
template<typename T, const int m, const int n>
void Matrix<T,m,n>::AssignFromList(AppendList *head){
    int i; int j;
    for (;head != nullptr;head=head->next){
        T temp = getValue(head->i,head->j);
        i = head->i;
        j = head->j;
        setValue(head->i,head->j, temp+(T)head->value);
    }

}

/// Display contents of the matrix
template<typename T, const int m, const int n>
void Matrix<T,m,n>::displayMatrix(){
    
    for (int i=0; i<m; i++){
        Vector<T,n> r;
        getRow(i,r);
        r.displayVector();
        std::cout << std::endl;
    }    

}

#endif