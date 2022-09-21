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

/* -------------------------------------------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------------------------------------------- */
template<typename T>
class Vector{
    private:
        T *m_Vector;
        T size;
    public:
        Vector(); // Constructor
        Vector(int Size); // Constructor
        Vector(T *v, int Size); // Constructor to initialize with values
        void init_(); // sets to null_ptr
        void setSize(int Size); // can explicitly set size
        int getLength_(){
            return size;
        }
        void setValue(int idx, T value); // sets value at some index
        T getValue(int idx); // gets value at that index
        void Add(T val); // adding scalar value to entire vector
        void Add(Vector<T> v_in, Vector<T> &v_out); // adding an entire vector to this vector and assign it to another vector
        void Subtract(T val, Vector<T> &v_out); 
        void Subtract(Vector<T> v_in, Vector<T> &v_out);
        void Scale(T scale, Vector<T> &v_out);
        void ElementMultiplication(Vector<T> v_in, Vector<T> &v_out);
        T dotProduct(Vector<T> &a);
        T *getData(int &Size);
        void AssignFromList(AppendList *head);
        void displayVector();
        //~Vector(){delete[] m_Vector;}
};

/// default constructor
template<typename T>
Vector<T>::Vector(){
//   init_();
}

/// @brief Constructor to initialize with its size and set its values to 0
/// @tparam T 
/// @param Size 
template<typename T>
Vector<T>::Vector(int Size){
    this->setSize(Size);
}

/// @brief Initialize new object vector to another pointer vector
/// @tparam T 
/// @param v - the pointer array whose values we set to the new Vector object
/// @param Size - Size of the array 
template<typename T>
Vector<T>::Vector(T *v, int Size){
    m_Vector = new T[Size];
    size = Size;
    for(int i=0; i<Size; i++){
        m_Vector[i] = v[i];
    }
}

/// @brief Initialize the object to a null_ptr
/// @tparam T 
template<typename T>
void Vector<T>::init_(){
    m_Vector = nullptr;
    size = 0;
}

/// @brief sets the size of the vector if not defined through constructor, 
/// maybe we can use this to reassign its size later if needed. Also assigns all values to 0
/// @tparam T 
/// @param Size size of the pointer vector
template<typename T>
void Vector<T>::setSize(int Size){
    m_Vector = nullptr;
    m_Vector = new T[Size];
    size = Size;
    for (int i=0; i<size; i++){
        m_Vector[i] = (T)(0.);
    }
}

/// @brief set value at index i 
/// @tparam T 
/// @param idx - index to set the value
/// @param value  - value we are going to assign
template<typename T>
void Vector<T>::setValue(int idx, T value){
    if(idx >= size){ std::cerr << "Index length greater than size" << std::endl;}
    else{ m_Vector[idx] = value; }
}

/// @brief Gets the value assigned to that index of this Vector
/// @tparam T 
/// @param idx - index for which we need the value assigned
/// @return 
template<typename T> 
T Vector<T>::getValue(int idx){
    if(idx>=size){std::cerr << "Index youre trying to access exceeds vector length" << std::endl; return (T)0.0;}
    else{ return m_Vector[idx]; }
}

template<typename T>
void Vector<T>::Add(T val){
    for(int i=0; i< size; i++){
        setValue(i,m_Vector[i]+val);
    }
}

/// @brief adds this vector to v_in and assigns it to v_out
/// @tparam T 
/// @param v_in - vector to be added to object
/// @param v_out - vector we need to assign the addition
template<typename T>
void Vector<T>::Add(Vector<T> v_in, Vector<T> &v_out){
    if(size == v_in.getLength_() && size == v_out.getLength_()){
        for (int i=0; i<size; i++){
            v_out.setValue(i, m_Vector[i]+v_in.getValue(i));
        }
    }
    else{
        std::cerr << "arrays of incompatible sizes cant be added \n";
    }
}

/// @brief subtracts val from each vector index and assigns it to v_out
/// @tparam T 
/// @param val - value to be subtracted FROM Vector object
/// @param v_out - Vector to be assigned to
template<typename T>
void Vector<T>::Subtract(T val, Vector<T> &v_out){
    if(size == v_out.getLength_()){
        for(int i=0; i<size; i++){
            v_out.setValue(i,m_Vector[i]-val);
        }
    }
    else{
        std::cerr << "incompatible vector lengths \n";
    }
}

/// @brief Subtract the current Vector object with the Vector v_in and assigns it to v_out, 
/// v_out = (current Vector object) - v_in
/// @tparam T 
/// @param v_in - vector to be subtracted FROM vector object
/// @param v_out - vector to be assigned this subtraction
template<typename T>
void Vector<T>::Subtract(Vector<T> v_in, Vector<T> &v_out){
    if (size == v_in.getLength_() && size == v_out.getLength_()){
        for(int i=0; i<size; i++){
            v_out.setValue(i,m_Vector[i]-v_in.getValue(i));
        }
    }
    else{
        std::cerr << "Incompatible Vector sizes \n";
    }
    
}

/// @brief Scales the current Vector object by (scale) and assigns it to Vector v_out
/// v_out = (scale)*(current Vector object)
/// @tparam T 
/// @param scale - scalar value used to scale
/// @param v_out - vector to assign after scaling Vector object
template<typename T>
void Vector<T>::Scale(T scale, Vector<T> &v_out){
    for(int i=0; i<size; i++){
        v_out.setValue(i,m_Vector[i] * scale);
    }
}

/// @brief Perform element by element multiplication, v_out[i]=(object)[i]*v_in[i]
/// @tparam T 
/// @param v_in -  vector used to multiply element vise
/// @param v_out - vector assigned the final value
template<typename T>
void Vector<T>::ElementMultiplication(Vector<T> v_in, Vector<T> &v_out){
    if (size == v_in.getLength_() && size==v_out.getLength_()){
        for(int i=0; i<size; i++){
            v_out.setValue(i, getValue(i)*v_in.getValue(i)); 
        }
    }
    else{
        std::cerr << "Incompatible vector sizes \n";
    }
    
}

/// @brief Performs v_in(.)a
/// @tparam T 
/// @param a - vector to dot product with
/// @return - returns scalar dot product value
template<typename T>
T Vector<T>::dotProduct(Vector<T> &a){
    T sum= (T)0.0;
    if (size == a.getLength_()){
        for (int i=0; i<size; i++){
            sum += a.getValue(i)*(this->getValue(i));
        }
        return sum;
    }
    else{
        std::cerr << "incompatible Vector sizes and Dot product cant be performed \n";
        return sum;
    }
    
}

/// gives the pointer to the place where m_Vector is stored. 
template<typename T>
T *Vector<T>::getData(int &Size){
    T *temp_vector = nullptr;
    temp_vector = m_Vector;
    Size = size;
    return temp_vector;
}

template<typename T>
void Vector<T>::AssignFromList(AppendList * head){
    int i;
    for (;head != nullptr;head=head->next){
        i = head->i;
        if (i < size){
            T temp = getValue(head->i);
            setValue(head->i, temp+(T)head->value);
        }
        else{
            std::cerr << "trying to reach non-existent index of this Vector \n";
        }
        
    }
}

/// Prints vector
template<typename T>
void Vector<T>::displayVector(){
    for(int i=0; i<size; i++){std::cout<<m_Vector[i]<<" ";}
}
/* --------------------------------------------------------------------------------------------------------------------- */
/* ---------------------------------------- Matrix class --------------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------------------------------- */

template<typename T>
class Matrix{
    private:
        T **m_Matrix;
        int m, n;
        bool SQUARE;
    public:
        Matrix();
        Matrix(const Matrix<T> *m_in);
        Matrix(int M, int N);
        Matrix(T *rowMajor, int Size, int M, int N); // imports matrix data from a rowMajor format
        void setSize(int M, int N);
        void getSize_(int &M, int &N){
            M = m; N = n;
        }
        Matrix(Vector<T> &v, Vector<T> &w);
        //Matrix(T i, int M); // Matrix to generate i*Identity matrix 
        void init_();
        void setValue(int i, int j, T val);
        T getValue(int i, int j);
        void setRow(int row_idx, Vector<T> &r);
        void setColumn(int col_idx, Vector<T> &c);
        void getRow(int row_idx, Vector<T> &r);
        void getColumn(int col_idx, Vector<T> &c);
        void Transpose(Matrix<T> &m_Matrix_T);
        void Add(T val, Matrix<T> &M_out);
        void Add(Matrix<T> M_in, Matrix<T> &M_out);
        void Multiply(T scale, Matrix<T> &M_out);
        void Multiply(Vector<T> u, Vector<T> &u_out);
        void Multiply(Matrix<T> M_in, Matrix<T> &M_out);
        void ElementMultiplication(Matrix<T> m1, Matrix<T> &m2);
        bool checkSquare(){return SQUARE;}
        void Invert(Matrix<T> &inv);
        T *exportRowMajor();
        void setRowMajor(T *rowMajor, int Size);
        void AssignFromList(AppendList *head);
        T **returnMatrix(){return m_Matrix;}
        void displayMatrix();
        // ~Matrix(){
        //     std::cout << "destructor array start \n";
        //     for (int i=0; i<m; i++){
        //         delete[] m_Matrix[i];
        //     }
        //     delete[] m_Matrix;
        //     std::cout << "destructor array end \n";
        // }
};

template<typename T>
Matrix<T>::Matrix(){
setSize(1,1);
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T> *m_in){
    m = m_in->m;
    n = m_in->n;
    SQUARE = m_in->SQUARE;
    m_Matrix = m_in->m_Matrix;
}

template<typename T>
Matrix<T>::Matrix(int M, int N){
    setSize(M,N);
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            m_Matrix[i][j] = (T)(0.);
        }
    }
}

/// Constructor to generate a Matrix object from a rowMajor array 
template<typename T>
Matrix<T>::Matrix(T *rowMajor, int Size, int M, int N){
    if (M*N != Size){
        std::cerr<<"Matrix size is incompatible" << std::endl;
    }
    else {
        setSize(M,N);
        int idx = 0;
        for (int i=0; i<M; i++){
            for (int j=0; j<N; j++){
                m_Matrix[i][j] = rowMajor[idx];
                ++idx;
                }
        }
    }
}

template<typename T>
void Matrix<T>::setSize(int M, int N){
    m_Matrix = nullptr;
    m_Matrix = new T *[M];
    for (int i=0; i< M; i++){
        m_Matrix[i] = new T[N];
    }
    m = M; n = N;
    SQUARE = (m==n)?true:false;
}

// /// Constructor to generate a Matrix object which is i*Identity matrix
// template<typename T>
// Matrix<T>::Matrix(T i, int M){
//     setSize(M,M);
//     for(int idx=0; idx<m; idx++){
//         m_Matrix[idx][idx] = (T)(i*1.);
//     }
// }


template<typename T>
Matrix<T>::Matrix(Vector<T> &v, Vector<T> &w){
    int M = v.getLength_();
    int N = w.getLength_();
    setSize(M,N);
    for (int i=0; i<m; i++){
        Vector<T> temp(N);
        T s = v.getValue(i);
        w.Scale(s,temp);
        setColumn(i,temp);
    }
}

/// Initialize all values with values 0.0
template<typename T>
void Matrix<T>::init_(){
    m_Matrix = new T *[1];
    for (int j=0; j<1; j++){
        m_Matrix[j] = nullptr;
    }
    m = 1; n = 0;
    SQUARE = true;
}

/// Set one single value of the current Matrix object
template<typename T>
void Matrix<T>::setValue(int i, int j, T val){
    if(i>=m || j>=n){
        std::cerr << "Index exceeds matrix dimensions" << std::endl;
    }
    else {
        m_Matrix[i][j] = val;
    }
}

/// Get value from a index i,j
template<typename T>
T Matrix<T>::getValue(int i, int j){
    if(i>=m || j>=n){
        std::cerr<<"Trying to access invalid index values" << std::endl;
        return (T)0.0;
    }
    else{
        return m_Matrix[i][j];
    }
}

/// Set row of current matrix of index row_idx with Vector r
template<typename T>
void Matrix<T>::setRow(int row_idx, Vector<T> &r){
    if (row_idx >= m){
        std::cerr << "Row index exceeds Matrix dimensions" << std::endl;
    }
    else if(r.getLength_() != n){
        std::cerr << "Vector length not equal to column length \n";
    }
    else {
        for(int i = 0; i < n; i++){
            m_Matrix[row_idx][i] = r.getValue(i);
        }
    }
}

/// Set column of current matrix of index col_idx with Vector c
template<typename T>
void Matrix<T>::setColumn(int col_idx, Vector<T> &c){
    if (col_idx >= n){
        std::cerr << "Column index exceeds Matrix dimensions" << std::endl;
    }
    else if(c.getLength_()!= m){
        std::cerr << "Vector length not equal to the row length \n";
    }
    else {
        for (int i = 0; i < m; i++){
            m_Matrix[i][col_idx] = c.getValue(i);
        }
    }
}

/// Get a row of current Matrix object and assign it to Vector r
template<typename T>
void Matrix<T>::getRow(int row_idx, Vector<T> &r){
    if(row_idx >= m){
        std::cerr << "Row index exceeds matrix dimensions" << std::endl;
    }
    else if(r.getLength_() != n){
        std::cerr << "Vector dimension not equal to column length \n";
    }
    else {
        for(int i=0; i<n; i++){
            r.setValue(i,m_Matrix[row_idx][i]);
        }
    }
}

/// Get a column of current Matrix object and assign it to Vector c
template<typename T>
void Matrix<T>::getColumn(int col_idx, Vector<T> &c){
    if(col_idx >= n){
        std::cerr << "Column index exceeds matrix dimensions" << std::endl;
    }
    else if(c.getLength_() != m){
        std::cerr << "Vector dimension not equal to row length \n";
    }
    else {
        for(int i=0; i<m; i++){
            c.setValue(i,m_Matrix[i][col_idx]);
        }
    }
}

/// Transposes current Matrix object and assigns it to m_Matrix_T
/// m_Matrix_T = (current Matrix)'
template<typename T>
void Matrix<T>::Transpose(Matrix<T> &m_Matrix_T){
    int mt,nt;
    m_Matrix_T.getSize_(mt,nt);
    if (mt == n && nt == m){
        for(int i=0; i<mt; i++){
            for(int j=0; j<nt; j++){
                m_Matrix_T.setValue(i,j,m_Matrix[j][i]);
            }
        }
    }
    else{
        std::cerr << "Matrix dimensions dont match up \n";
    }
    
}

/// Add a scalar value (val) and add it to all elements of current Matrix object
/// and assigns it to M_out. M_out = (current Matrix) + val
template<typename T>
void Matrix<T>::Add(T val, Matrix<T> &M_out){
    int mo, no;
    M_out.getSize_(mo,no);
    if (mo == m && no == n){
        for(int i=0; i<m; i++){
            Vector<T> r(n);
            getRow(i,r);
            r.Add(val,r);
            M_out.setRow(i,r);
        }
    }
    else{
        std::cerr << "Matrix dimensions dont match \n";
    }
    
}

/// Adds current matrix object to Matrix object M_in and assigns to M_out
/// M_out = (current Matrix) + M_in
template<typename T>
void Matrix<T>::Add(Matrix<T> M_in, Matrix<T> &M_out){
int mi, ni, mo, no;
M_in.getSize_(mi,ni); 
M_out.getSize_(mo,no);
if (mi==m && ni==n && mo==m && no==n){
    T *rm1 = exportRowMajor();
    T *rm2 = M_in.exportRowMajor();
    T *vec = (T*)malloc((m*n)*sizeof(T));
    for (int i=0; i < m*n; i++){
        vec[i] = rm1[i] + rm2[i];
    }
    Matrix<T> test(vec, m*n, m, n);
    //M_out.setRowMajor(vec, m*n);
    M_out = test;    
}
else{
    std::cerr << "Matrix dimensions dont match \n";
}
}

/// Scales each element of current Matrix object by (value scale) and assigns it to
/// Matrix object M_out. M_out = (scale)*(currentMatrix)
template<typename T>
void Matrix<T>::Multiply(T scale, Matrix<T> &M_out){
    int mo, no;
    M_out.getSize_(mo,no);
    if (m==mo && n==no){
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                M_out.setValue(i,j,m_Matrix[i][j]*scale);
            }
        }        
    }
    else{
        std::cerr << "Matrix dimensions are incompatible \n";
    }
}

/// This multiplies current Matrix object with a vector u and assigns it to 
/// another vector u_out. u_out = (currentMatrix)*u 
template<typename T>
void Matrix<T>::Multiply(Vector<T> u, Vector<T> &u_out){
    int ni = u.getLength_();
    int mo = u_out.getLength_();

    if (n == ni && mo == m){
        for(int i=0; i<m; i++){
            Vector<T> row(n) ;
            getRow(i,row);
            T val = u.dotProduct(row);
            u_out.setValue(i,val);
        }
    }
    else{
        std::cerr << "Matrix dimensions are incompatible \n";
    }

}

/// Current Matrix object is matrix multiplied with Matrix object M_in and assigned to
/// Matrix object M_out. M_out = (current)*M_in 
template<typename T>
void Matrix<T>::Multiply(Matrix<T> M_in, Matrix<T> &M_out){
int mi, ni;
M_in.getSize_(mi,ni);
if (n == mi){
    for(int i=0; i<m; i++){
        Vector<T> row(n);
        getRow(i,row);
        for(int j=0; j<ni; j++){
            Vector<T> column(ni);
            M_in.getColumn(j,column);
            T val = row.dotProduct(column);
            M_out.setValue(i,j,val);
        }
    }    
}
else{
    std::cerr << "Matrix dimensions are incompatible \n";
}
}

/// Element by Element multiplication is performed on current Matrix object and Matrix
/// object M_in and assigns it to Matrix object M_out
template<typename T>
void Matrix<T>::ElementMultiplication(Matrix<T> M_in, Matrix<T> &M_out){
    Vector<T> temp1(m), temp2(m), temp3(m);
    for(int j=0; j<n; j++){
        getColumn(j,temp1);
        M_in.getColumn(j,temp2);
        temp1.ElementMultiplication(temp2,temp3);
        M_out.setColumn(j,temp3);
    }
}

/// Inverts the current matrix and assigns it to Matrix object inverse
template<typename T>
void Matrix<T>::Invert(Matrix<T> &inv){
    if(SQUARE){

    }
    else {
        std::cerr<< "Matrix is not a square matrix and cannot be inverted" << std::endl;
    }
}

/// exports the matrix as a pointer to an array. The rows of the matrix are organized one next to each other.
template<typename T>
T *Matrix<T>::exportRowMajor(){
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
template<typename T>
void Matrix<T>::setRowMajor(T *rowMajor, int Size){
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
template<typename T>
void Matrix<T>::AssignFromList(AppendList *head){
    int i; int j;
    for (;head != nullptr;head=head->next){
        T temp = getValue(head->i,head->j);
        i = head->i;
        j = head->j;
        setValue(head->i,head->j, temp+(T)head->value);
    }

}

/// Display contents of the matrix
template<typename T>
void Matrix<T>::displayMatrix(){
    
    for (int i=0; i<m; i++){
        Vector<T> r(n);
        getRow(i,r);
        r.displayVector();
        std::cout << std::endl;
    }    

}

#endif