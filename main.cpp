#include <iostream>
#include <fstream>

//using namespace cimg_library;

class Point {
public:
    Point(){
        properties_count = 0;
    }
    Point(int* x, int count){
        std::cout << "point constructor\n";
        this->properties_count = count;
        this->properties = new int[count];
        for (int i = 0; i < count; i++){
            properties[i] = x[i];
        }
    }
    int get_properties(int i){
        return properties[i];
    }
    int get_size(){
        return properties_count;
    }
    int* get_all(){
        return properties;
    }
private:
    int* properties;
    int properties_count;
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
    Matrix(int rows_count, Point* points){
        std::cout << "matrix constructor\n";
        n = rows_count;
        m = 0;
        for (int i = 0; i < n; i++){
            m = max(m, points[i].get_size());
        }
        values = new double*[n];
        for (int i = 0; i < n; i++){
            values[i] = new double[m];
            for (int j = 0; j < m; j++){
                if (j < points[i].get_size()){
                    values[i][j] = points[i].get_properties(j);
                }else{
                    values[i][j] = 0;
                }
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
                std::cout << values[i][j];
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
    SquareMatrix(int rows_count, Point *points) : Matrix(rows_count, points) {
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

Matrix operator * (Matrix a, int b){
    Matrix new_mat(a);
    for (int i = 0; i < a.get_size_n(); i++){
        for (int j = 0; j < a.get_size_m(); j++){
            new_mat.get_values()[i][j] *= b;
        }
    }
    return new_mat;
}

Matrix operator * (Matrix a, double b){
    Matrix new_mat(a);
    for (int i = 0; i < a.get_size_n(); i++){
        for (int j = 0; j < a.get_size_m(); j++){
            new_mat.get_values()[i][j] *= b;
        }
    }
    return new_mat;
}

SquareMatrix operator * (SquareMatrix a, int b){
    SquareMatrix new_mat(a);
    for (int i = 0; i < a.get_size_n(); i++){
        for (int j = 0; j < a.get_size_m(); j++){
            new_mat.get_values()[i][j] *= b;
        }
    }
    return new_mat;
}

SquareMatrix operator * (SquareMatrix a, double b){
    SquareMatrix new_mat(a);
    for (int i = 0; i < a.get_size_n(); i++){
        for (int j = 0; j < a.get_size_m(); j++){
            new_mat.get_values()[i][j] *= b;
        }
    }
    return new_mat;
}

Matrix operator * (Matrix a, Matrix b){
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
    }
}
Matrix linearReg(Matrix x, Matrix y){
    Matrix m_1 = (x.get_trans_matrix() * x);
    m_1.show();
    SquareMatrix m_2 = SquareMatrix(m_1);

    std::cout << "det = " << m_2.get_det() << "\n";
    if (m_2.get_det() != 0){
        Matrix m_3 = m_2.get_reverse_matrix() * x.get_trans_matrix();
        m_3.show();
        m_3 = m_3 * y.get_trans_matrix();
        m_3.show();
        return  m_3;
    }
    return m_2;
}
/*
Point** readMyPoints(FILE* f){
    std::ifstream file;
    file.open("../files/f_1.txt");
    if (file.is_open()){
        Point** ans =
        int n, m;
        file >> n >> m;
        for (int i = 0; i < n; i++){
            int p_properties[m];
            for (in j = 0; j < m; j++){
                file >> p_properties[j];
            }
            Point p(p_properties, m);
        }
    }
    file.close()
    }
}//*/
/*
void writeToFPoints(){
    std::ofstream file;
    file.open("../files/f_1.txt");
    if (file.is_open()){
        int points_count = 10, properties_count = 2;
        int v_x[] = {10, 5, 12, 25, 1, 18, 11, 7, 19, 15};
        file << points_count << " " << properties_count << "\n";
        for (int i = 0; i < points_count - 1; i++) {
            file << i << " " << v_x[i] << ",";
        }
        file << (points_count - 1) << " " << v_x[points_count - 1] << "\n";

        int v_y[] = {15, 12, 18, 30, 3, 20, 14, 10, 20, 13};
        for (int i = 0; i < points_count - 1; i++) {
            file << i << " " << v_y[i] << ",";
        }
        file << (points_count - 1) << " " << v_y[points_count - 1];
    }
    file.close();
}//*/

int main() {
    std::cout << "Hello, World!" << std::endl;
    int v_x[] = {10, 5, 12, 25, 1, 18, 11, 7, 19, 15},
        v_y[] = {15, 12, 18, 30, 3, 20, 14, 10, 20, 13};
    Point p_y(v_y, 10);
    Point px[] = {Point(v_x, 10)};
    Point py[] = {p_y};
    Matrix x(1, px), y(1, py);
    x.show();
    y.show();
    //linearReg(x, y);

    writeToFPoints();

    return 0;
}
