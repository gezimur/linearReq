#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

//using namespace cimg_library;

class Point {
public:
    Point(){
        properties_count = 0;
    }
    Point(const double* x, int count, double y){
        std::cout << "point constructor\n";
        this->properties_count = count;
        this->properties = new double[count];
        for (int i = 0; i < count; i++){
            properties[i] = x[i];
        }
        this->y = y;
    }

    double* properties;
    int properties_count;
    double y;
};

class Matrix{
public:
    Matrix(int n, int m){
        this->n = n;
        this->m = m;
        values = new double*[n];
        for (int i = 0; i < n; i++){
            values[i] = new double[m];
        }
    }
    Matrix(Point* points, int points_count, bool is_x_matrix){
        n = points_count;
        values = new double*[n];
        if (is_x_matrix){
            m = points[0].properties_count;
            for (int i = 0; i < n; i++){
                values[i] = new double[m];
                for (int j = 0; j < m; j++){
                    values[i][j] = points[i].properties[j];
                }
            }
        }else{
            m = 1;
            for (int i = 0; i < n; i++){
                values[i] = new double[m];
                values[i][0] = points[i].y;
            }
        }
    }
    Matrix(Matrix const &mat){
        n = mat.n;
        m = mat.m;
        values = new double*[n];
        for (int i = 0; i < n; i++){
            values[i] = new double[m];
            for (int j = 0; j < m; j++){
                values[i][j] = mat.values[i][j];
            }
        }
    }

    double** get_values(){
        return values;
    }
    double get_value(int i, int j){
        return values[i][j];
    }
    int get_size_n(){
        return n;
    }
    int get_size_m(){
        return m;
    }
    double get_average_in_row(int i){
        double ans = values[i][0];
        for (int j = 1; j < m; j++){
            ans += values[i][j];
        }
        return ans / m;
    }

    Matrix get_trans_matrix(){
        Matrix new_mat(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                new_mat.get_values()[i][j] = values[j][i];
            }
        }
        return new_mat;
    }

    void show(){
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                std::cout.width(6);
                std::cout.precision(3);
                std::cout << values[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

protected:
    int n,m;
    double** values;

    static int max(int a, int b) {
        return (a > b) ? a : b;
    }

};

class SquareMatrix :public Matrix{
public:
    SquareMatrix(int n, int m) : Matrix(n, m){
        std::cout << "SquareMatrix constructor\n";
    }

    SquareMatrix(Point *points, int n, bool is_x_matrix) : Matrix(points, n, is_x_matrix) {
        std::cout << "SquareMatrix constructor\n";
    }

    explicit SquareMatrix(Matrix& mat) : Matrix(mat) {
        std::cout << "SquareMatrix constructor\n";
    }

    double get_det(){
        double ans = 0;
        bool columns[m];
        for (int i = 0; i < m; i++){
            columns[i] = true;
        }
        for (int i = 0; i < m; i++){
            if (i % 2 == 0){
                ans += _get_A(0, i, columns, n + 1,m + 1);
            }else{
                ans -= _get_A(0, i, columns, n + 1, m + 1);
            }
        }
        return ans;
    }

    double get_A(int a_i, int a_j){
        double ans = 0;
        bool columns[m];
        for (int i = 0; i < m; i++){
            columns[i] = true;
        }
        columns[a_j] = false;
        int z = 1;
        for (int i = 0; i < m; i++){
            if (columns[i]){
                ans += z * _get_A(0, i, columns, a_i,a_j);
                z *= -1;
            }
        }
        z = ((a_i + a_j) % 2 == 0)? 1 : -1;
        return z * ans;
    }

    SquareMatrix get_reverse_matrix(){
        double det = get_det();
        if (det != 0){
            SquareMatrix new_mat(n, m);
            for (int i = 0; i < n; i++){
                for (int j = 0; j < m; j++){
                    new_mat.get_values()[i][j] = get_A(i, j);
                }
            }
            //new_mat.show();
            return (new_mat * (1.0 / det)).get_trans_squarematrix();
        }
    }

    SquareMatrix get_trans_squarematrix(){
        SquareMatrix new_mat(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                new_mat.get_values()[i][j] = values[j][i];
            }
        }
        return new_mat;
    }

    SquareMatrix operator * (double b){
        SquareMatrix new_mat(*this);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                new_mat.get_values()[i][j] *= b;
            }
        }
        return new_mat;
    }

private:
    double _get_A(int i, int j, bool* columns, int a_i, int a_j){
        if (i == a_i) i++;
        if (i == n - 1){
            return values[i][j];
        }/*
        if (i > n - 1){
            return 1;
        }//*/
        columns[j] = false;
        double ans = 0;
        int z = 1;
        bool flag = true;
        for (int t = 0; t < m; t++){
            if (columns[t] && t != a_j){
                ans += z * values[i][j] * _get_A(i + 1, t, columns, a_i ,a_j);
                z *= -1;
                flag = false;
            }
        }
        if (flag) ans = values[i][j];
        columns[j] = true;
        return ans;
    }
};

Matrix operator * (Matrix &a, int &b){
    Matrix new_mat(a);
    for (int i = 0; i < a.get_size_n(); i++){
        for (int j = 0; j < a.get_size_m(); j++){
            new_mat.get_values()[i][j] *= b;
        }
    }
    return new_mat;
}

Matrix operator * (Matrix &a, double &b){
    Matrix new_mat(a);
    for (int i = 0; i < a.get_size_n(); i++){
        for (int j = 0; j < a.get_size_m(); j++){
            new_mat.get_values()[i][j] *= b;
        }
    }
    return new_mat;
}

SquareMatrix operator * (SquareMatrix &a, int &b){
    SquareMatrix new_mat(a);
    for (int i = 0; i < a.get_size_n(); i++){
        for (int j = 0; j < a.get_size_m(); j++){
            new_mat.get_values()[i][j] *= b;
        }
    }
    return new_mat;
}

SquareMatrix operator * (SquareMatrix &a, double &b){
    SquareMatrix new_mat(a);
    for (int i = 0; i < a.get_size_n(); i++){
        for (int j = 0; j < a.get_size_m(); j++){
            new_mat.get_values()[i][j] *= b;
        }
    }
    return new_mat;
}

Matrix operator * (Matrix &a, Matrix &b){
    if (a.get_size_m() == b.get_size_n()){
        int new_mat_n = a.get_size_n(), new_mat_m = b.get_size_m();
        SquareMatrix new_mat(new_mat_n, new_mat_m);
        for (int i = 0; i < new_mat_n; i++){
            for (int j = 0; j < new_mat_m; j++){
                double ans = 0;
                for (int t = 0; t < a.get_size_m(); t++){
                    ans += a.get_values()[i][t] * b.get_values()[t][j];
                }
                new_mat.get_values()[i][j] = ans;
            }
        }
        return new_mat;
    }else{
        std::cout << "a.get_size_m() != b.get_size_n() \n";
        return Matrix(0,0);
    }
}

int linearReg(Matrix &x, Matrix &y){
    Matrix m = x.get_trans_matrix();
    //std::cout << "xt :\n";
    //m.show();

    m = (m * x);
    //std::cout << "xt * x :\n";
    //m.show();

    SquareMatrix sm = SquareMatrix(m).get_reverse_matrix();
    //std::cout << "(xt * x)^(-1) :\n";
    //sm.show();

    if (sm.get_det() != 0){
        //std::cout << "det(xt * x) =" << sm.get_det() << "\n";
        m = x.get_trans_matrix();

        Matrix m_add = sm * m;
        //std::cout << "(xt * x)^(-1) * xt :\n";
        //m_add.show();

        m = y;
        //std::cout << "y :\n";
        //m.show();

        m_add = m_add * m;
        //std::cout << "O :\n";
        m_add.show();
        return 0;
    }else{
        std::cout << "det = 0\n";
        return 1;
    }
}

Point* createTestPoints(int n, int m){
    srand( time(nullptr) );
    Point* ans = new Point[n];
    for (int i = 0; i < n; i++){
        double x[m + 1];
        double s = 0;
        for (int j = 1; j < m + 1; j++) {
            x[j] = rand() % 100 + (rand() % 100) / 100.0;
            s += x[j];
        }
        x[0] = 1;
        double y = sin(s) * 5 + s * 0.3 + 5;
        ans[i] = Point(x, m + 1, y);
    }
    return ans;
}

void wirtePoints(Point* points, int n){
    std::ofstream f;
    f.open("../files/f_1.txt");
    for (int i = 0; i  < n; i++){
        for (int j = 1; j < points[i].properties_count; j++){
            f << points[i].properties[j] << " ";
        }
        f << points[i].y << "\n";
    }
    f.close();
}

int main() {
    std::cout << "Hello, World!" << std::endl;

    //linearReg(x, y);
    int n = 10, m = 1;
    Point* points = createTestPoints(n, m);
    Matrix x(points, n, true), y(points, n, false);

    //std::cout << "x :\n";
    //x.show();
    //std::cout << "y :\n";
    //y.show();
    linearReg(x, y);

    wirtePoints(points, n);

    return 0;
}
