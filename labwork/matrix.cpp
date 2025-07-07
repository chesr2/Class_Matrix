#include "matrix.h"
namespace linalg {
	const double eps = 0.000000001;

	//Конструкторы
	
	//Определение перегрузки конструктора Matrix, obj(size_t ,size_t )*
	Matrix::Matrix(int rows, int columns)
		:m_rows(rows), m_columns(columns)
	{
		if (rows<0||columns<0)
			throw std::runtime_error("Error, value of rows or columns cannot be minus");

		if (m_rows != 0 || m_columns != 0) {
			m_ptr = new double[m_rows * m_columns];
			for (size_t i = 0; i < m_rows * m_columns; ++i)
				m_ptr[i] = 0;
		}
		else
			m_ptr = nullptr;
	}

	Matrix::Matrix(int rows) {
		if (rows != 0) {
			m_rows = rows;
			m_columns = 1;
			m_ptr = new double[m_columns*m_rows];
			for (size_t i = 0; i < m_rows * m_columns; ++i)
				m_ptr[i] = 0;
		}
	}

	//Определение перегрузки конструктора Matrix, Matrix(const Matrix& )\
	(присваивающее копирование)*
	Matrix::Matrix(const Matrix& obj)
		:Matrix(obj.m_rows, obj.m_columns)
	{
		if (m_ptr != nullptr)
			for (size_t i = 0; i < m_rows * m_columns; ++i)
				m_ptr[i] = obj.m_ptr[i];
		else
			m_ptr = nullptr;
	}

	//Определение перегрузки конструктора Matrix, Matrix(Matrix&& )\
	(перемещающее копирование)*
	Matrix::Matrix(Matrix&& obj) noexcept {
		std::swap(m_ptr, obj.m_ptr);
		std::swap(m_rows, obj.m_rows);
		std::swap(m_columns, obj.m_columns);
	}

	//Опредение перегрузки конструктора Matrix, obj = {...}*
	Matrix::Matrix(std::initializer_list<double> lst)
		:Matrix(lst.size())
	{
		size_t i = 0;
		for (double el : lst)
			m_ptr[i++] = el;
	}

	//Определение перегрузки конструктора Matrix, obj = {{...},{...}...}*
	Matrix::Matrix(std::initializer_list<std::initializer_list<double>> lst)
		:Matrix(lst.size(), lst.begin()->size())
	{
		bool Usl = true;
		size_t h = lst.begin()->size();
		for (std::initializer_list<double> el_lst : lst)
		{
			if (h != el_lst.size())
				Usl = false;
		}
		size_t i = 0;
		if (Usl)
			for (std::initializer_list<double> el_lst : lst)
			{
				for (double el_lst_lst : el_lst)
				{
					m_ptr[i++] = el_lst_lst;
				}
			}
		else
			throw std::runtime_error("Error, sizes of rows in this list are not same");
	}

	//Функции-методы

	//Определение функции-метода .empty()*
	const bool Matrix::empty() const noexcept {
		if (m_ptr == nullptr && m_rows == 0 && m_columns == 0)
			return true;
		return false;
	}

	//Определение функции-метода .reshape()*
	void Matrix::reshape(const size_t rows, const size_t columns) {
		if (rows * columns != m_columns * m_rows)
			throw std::runtime_error("Error, cannot reshape because of incorrect sizes, elements are deleting");

		m_rows = rows;
		m_columns = columns;
	}

	//Определение функции-метода вычисляющей норму Фробениуса*
	const double Matrix::norm() const {
		double norm = 0;
		if (m_ptr != nullptr) {
			for (size_t i = 0; i < m_columns * m_rows; ++i) {
				norm += m_ptr[i] * m_ptr[i];
			}
			return sqrt(norm);
		}
		else
			throw std::runtime_error("Error, matrix is empty");
		return 0;
	}

	//Определение функции-метода вычисляющей след*
	const double Matrix::trace() const {
		if (m_columns != m_rows)
			throw std::runtime_error("Error, not square matrix");
		double trace = 0;
		if (m_ptr != nullptr) {
			for (size_t i = 0; i < m_columns * m_rows; ++i) {
				if (i % m_columns == i / m_columns)
					trace += m_ptr[i];
			}
			return trace;
		}
		else
			throw std::runtime_error("Error, matrix is empty");
		return 0;
	}

	//Определение функции-метода вычисляющей определитель*
	const double Matrix::det() const{
		if (m_columns != m_rows)
			throw std::runtime_error("Error, not square matrix");

		Matrix obj=(*this);
		size_t i = 0;
		size_t value_swap = 0;

		for (size_t j = 0; i != m_rows - 1; ++j) {
			// Поиск главного элемента
			size_t maxRowIndex = i;
			for (size_t k = i + 1; k < m_rows; ++k) {
				if (std::abs(obj.m_ptr[k * m_columns + j]) > std::abs(obj.m_ptr[maxRowIndex * m_columns + j])) {
					maxRowIndex = k;
				}
			}

			// Перестановка строк, если необходимо
			if (maxRowIndex != i) {
				(obj).swap_rows(i, maxRowIndex);
				++value_swap;
			}

			// Приведение матрицы к ступенчатому виду
			if (obj.m_ptr[i * m_columns + j] != 0) {

				for (size_t k = i + 1; k < m_rows; ++k) {
					double temp = obj.m_ptr[k * m_columns + j] / obj.m_ptr[i * m_columns + j];

					for (size_t h = j; h < m_columns; ++h) {
						obj.m_ptr[k * m_columns + h] += (-1) * obj.m_ptr[i * m_columns + h] * temp;
						if (std::fabs(obj.m_ptr[k * m_columns + h]) < eps)
							obj.m_ptr[k * m_columns + h] = 0;
					}
				}

			}
			else {
				for (size_t k = i + 1; k < m_rows; ++k) {
					obj.m_ptr[k * m_columns + j] = 0;
				}
			}
			++i;
		}

		double det = 1;
		for (size_t i = 0; i < m_columns * m_rows; ++i)
			if (i / m_columns == i % m_columns) {
				det *= obj.m_ptr[i];
			}
		if(value_swap%2==0)
			return det;
		else
			return (-1)*det;
	}

	//Определение функции-метода находящей ранг матрицы*
	const size_t Matrix::rank() const{
		bool bl = 1;
		size_t rank = 0;
		Matrix obj = *this;
		obj.gauss_forward();
		for (size_t i = 0; i < obj.m_rows; ++i) {
			for (size_t j = 0; j < obj.m_columns; ++j) {
				if (std::fabs(obj.m_ptr[i * m_columns + j]) > eps) {
					bl = 0; break;
				}
			}
			if (!bl)
				++rank;
			bl = 1;
		}
		return rank;
	}

	//Определение функции-метода реализующей прямой метод Гаусса*
	Matrix& Matrix::gauss_forward() {
		size_t i = 0;
		size_t value_swap=0;

		for (size_t j = 0; j<m_columns&&i<m_rows; ++j) {
			// Поиск главного элемента
			size_t maxRowIndex = i;
			for (size_t k = i + 1; k < m_rows; ++k) {
				if (std::abs(m_ptr[k * m_columns + j]) > std::abs(m_ptr[maxRowIndex * m_columns + j])) {
					maxRowIndex = k;
				}
			}

			// Перестановка строк, если необходимо
			if (maxRowIndex != i) {
				(*this).swap_rows(i, maxRowIndex);
				++value_swap;
			}

			// Приведение матрицы к ступенчатому виду
			if (m_ptr[i * m_columns + j] != 0) {

				for (size_t k = i + 1; k < m_rows; ++k) {
					double temp = m_ptr[k * m_columns + j] / m_ptr[i * m_columns + j];

					for (size_t h = j; h < m_columns; ++h) {
						m_ptr[k * m_columns + h] += (-1) * m_ptr[i * m_columns + h] * temp;
						if (std::fabs(m_ptr[k * m_columns + h]) < eps)
							m_ptr[k * m_columns + h] = 0;
					}
				}

			}
			else {
				for (size_t k = i + 1; k < m_rows; ++k) {
					m_ptr[k * m_columns + j] = 0;
				}
			}
			++i;
		}
		return *this;
	}

	//Определение функции-метода реализующей обратный метод Гаусса*
	Matrix& Matrix::gauss_backward() {
		int iter = m_rows - 1;
		double temp = 0;

		//Умножаем строку на temp и вычитаем ее из следующей, цикл убывающий, потом работаем с матрицей i-1 j-1 и тд.
		for (int i = m_rows - 1; i > 0; --i) {
			for (int k = i - 1; k >= 0; --k) {
				if (m_ptr[i * m_columns + iter] != 0)
					temp = m_ptr[k * m_columns + iter] / m_ptr[i * m_columns + iter];
				else
					temp = 0;
				for (int h = m_columns - 1; h >= 0; --h) {
					m_ptr[k * m_columns + h] += (-1) * (m_ptr[i * m_columns + h]) * temp;
					if (std::fabs(m_ptr[k * m_columns + h]) < eps)
						m_ptr[k * m_columns + h] = 0;
				}
			}
			--iter;
		}
		return *this;
	}

	//Определение функции-метода удаляющей строку из матрицы*
	Matrix& Matrix::del_row(const size_t row) {
		Matrix obj(m_rows - 1, m_columns);
		size_t iter = 0;

		for (size_t i = 0; i < m_rows * m_columns; ++i) {
			if (i / m_columns != row) {
				obj.m_ptr[iter] = m_ptr[i];
				++iter;
			}
		}

		(*this) = obj;
		return *this;
	}

	//Определение функции-метода удаляющей столбец из матрицы*
	Matrix& Matrix::del_column(const size_t column) {
		Matrix obj(m_rows, m_columns - 1);
		size_t iter = 0;

		for (size_t i = 0; i < m_rows * m_columns; ++i) {
			if (i % m_columns != column) {
				obj.m_ptr[iter] = m_ptr[i];
				++iter;
			}
		}

		(*this) = obj;
		return *this;
	}

	//Определение функции-метода находящей минор элемента матрицы*
	const double Matrix::minor_ext(const size_t row, const size_t column) {
		Matrix obj = *this;
		obj.del_row(row);
		obj.del_column(column);
		return obj.det();
	}

	//Определение функции-метода строющей матрицу миноров*
	Matrix Matrix::mat_minor_ext() {
		Matrix obj(m_rows, m_columns);

		for (size_t i = 0; i < m_rows * m_columns; ++i) {
			obj.m_ptr[i] = (*this).minor_ext(i / m_columns, i % m_columns);
		}
		return obj;
	}

	//Определение функции-метода строющей матрицу алгебраических дополнений*
	Matrix Matrix::mat_alg_ext() {
		Matrix obj = (*this).mat_minor_ext();

		for (size_t i = 0; i < m_rows; ++i) {
			for (size_t j = 0; j < m_columns; ++j) {
				if ((i + j) % 2 != 0)
					obj.m_ptr[i * m_columns + j] *= (-1);
			}
		}
		return obj;
	}

	//Определение функции-метода меняющей местами строки матрицы*
	Matrix& Matrix::swap_rows(const size_t row1, const size_t row2) {
		double temp;

		for (size_t i = 0; i < m_columns; ++i) {
			temp = m_ptr[row1 * m_columns + i];
			m_ptr[row1 * m_columns + i] = m_ptr[row2 * m_columns + i];
			m_ptr[row2 * m_columns + i] = temp;
		}

		return *this;
	}

	//Определение функции-метода которая складывает все элементы матрицы*
	const double Matrix::sum_el() {
		double sum = 0;

		for (size_t i = 0; i < m_columns * m_rows; ++i) {
			sum += m_ptr[i];
		}
		return sum;
	}

	//Операторы

	//Определение перегрузки оператора = (присваивающее копирование)*
	Matrix& Matrix::operator = (const Matrix& obj) {
		m_rows = obj.m_rows;
		m_columns = obj.m_columns;

		if (obj.m_ptr != nullptr) {
			m_ptr = new double[m_rows * m_columns];
			for (size_t i = 0; i < m_rows * m_columns; ++i)
				m_ptr[i] = obj.m_ptr[i];
		}
		else
			m_ptr = nullptr;
		return *this;
	}

	//Определение перегрузки оператора = (перемещающее копирование)*
	Matrix& Matrix::operator = (Matrix&& obj) noexcept {
		std::swap(m_rows, obj.m_rows);
		std::swap(m_columns, obj.m_columns);
		std::swap(m_ptr, obj.m_ptr);
		return *this;
	}

	//Определение перегрузки оператора ()*
	const double Matrix::operator () (const size_t row, const size_t column) const {
		if (row >= m_rows || column >= m_columns) {
			throw std::runtime_error("Invalid indexes");
		}

		return m_ptr[column + m_columns * row];
	}

	double& Matrix::operator () (const size_t row, const size_t column) {
		if (row >= m_rows || column >= m_columns) {
			throw std::runtime_error("Invalid indexes");
		}

		return m_ptr[column + m_columns * row];
	}

	//Определение перегрузки оператора += (Поэлементное сложение МАТРИЦ)*
	Matrix& Matrix::operator += (const Matrix& obj) {
		if (m_rows != obj.m_rows || m_columns != obj.m_columns)
			throw std::runtime_error("Error, Size of matrixes are not same");

		if (m_ptr != nullptr) {
			for (size_t i = 0; i < m_columns * m_rows; ++i) {
				m_ptr[i] += obj.m_ptr[i];
			}
		}
		return *this;
	}

	//Определение перегрузки оператора -= (Поэлементное вычитание МАТРИЦ)*
	Matrix& Matrix::operator -= (const Matrix& obj) {
		if (m_rows != obj.m_rows || m_columns != obj.m_columns)
			throw std::runtime_error("Error, Size of matrixes are not same");

		if (m_ptr != nullptr) {
			for (size_t i = 0; i < m_columns * m_rows; ++i) {
				m_ptr[i] -= obj.m_ptr[i];
			}
		}
		return *this;
	}

	//Определение перегрузки оператора *= (перемножение МАТРИЦ)*
	Matrix& Matrix::operator *= (const Matrix& obj) {
		if (m_columns != obj.m_rows)
			throw std::runtime_error("This matrixes mustnt multiplicate");

		linalg::Matrix obj1(m_rows, obj.m_columns), obj_row, obj_column, objs_mlt;
		for (size_t i = 0; i < obj1.m_rows * obj1.m_columns; ++i) {
			obj_row = get_row(i / obj1.m_columns, *this);
			obj_row = transpose(obj_row);
			obj_column = get_column(i % obj1.m_columns, obj);
			for (size_t j = 0; j < obj_column.m_columns * obj_column.m_rows; ++j) {
				obj_column.m_ptr[j] *= obj_row.m_ptr[j];
			}
			obj1(i / obj1.m_columns, i % obj1.m_columns) = (obj_column).sum_el();
			if (std::fabs(obj1(i / obj1.m_columns, i % obj1.m_columns)) < eps)
				obj1(i / obj1.m_columns, i % obj1.m_columns) = 0;
		}
		*this = obj1;
		return *this;
	}

	//Определение перегрузки оператора * (Поэлементное перемножение на ЧИСЛО)*
	Matrix& Matrix::operator *= (const double val) {
		if (m_ptr != nullptr) {
			for (size_t i = 0; i < m_columns * m_rows; ++i) {
				m_ptr[i] *= val;
			}
		}
		return *this;
	}

	//Определение перегрузки оператора + (Поэлементное сложение МАТРИЦ)*
	const Matrix Matrix::operator + (const Matrix& obj) const {
		Matrix temp(*this);
		return temp += obj;
	}

	//Определение перегрузки оператора - (Поэлементное вычитание МАТРИЦ)*
	const Matrix Matrix::operator - (const Matrix& obj) const {
		Matrix temp(*this);
		return temp -= obj;
	}

	//Определение перегрузки оператора * (перемножение МАТРИЦ)*
	const Matrix Matrix::operator * (const Matrix& obj) const {
		Matrix temp(*this);
		return temp *= obj;
	}
	//Определение перегрузки оператора * (Поэлементное перемножение на ЧИСЛО)*
	const Matrix Matrix::operator * (const double val) const {
		Matrix temp(*this);
		return temp *= val;
	}

	//Определение перегрузки оператора * (Поэлементное перемножение на число)*
	const linalg::Matrix operator * (const double val, const linalg::Matrix& obj) {
		linalg::Matrix temp(obj);
		return temp *= val;
	}

	//Определение перегрузки оператора == (Проверка матриц на совпадение)*
	const bool operator == (const Matrix& obj1, const Matrix& obj2) {
		if (obj1.rows() != obj2.rows() || obj1.columns() != obj2.columns())
			return false;
		else {
			size_t i_row = 0, i_column = 0;
			for (size_t i = 0; i < obj1.columns() * obj1.columns(); ++i) {
				i_row = i / obj1.columns();
				i_column = i % obj1.columns();
				if (std::fabs(obj1(i_row, i_column) - obj2(i_row, i_column))>eps) {
					return false;
				}
			}
			return true;
		}
	}

	//Определение перегрузки оператора != (Проверка матриц на не совпадение)*
	const bool operator != (const Matrix& obj1, const Matrix& obj2) {
		if (obj1.rows() != obj2.rows() || obj1.columns() != obj2.columns())
			return true;
		else {
			size_t i_row = 0, i_column = 0;
			for (size_t i = 0; i < obj1.columns() * obj1.columns(); ++i) {
				i_row = i / obj1.columns();
				i_column = i % obj1.columns();
				if (std::fabs(obj1(i_row, i_column) - obj2(i_row, i_column)) > eps) {
					return true;
				}
			}
			return false;
		}
	}

	//Функция транспонирования*
	const Matrix transpose(const Matrix obj) {
		Matrix obj1 = obj;
		obj1.reshape(obj1.columns(), obj1.rows());
		for (size_t j = 0; j < obj1.rows(); ++j) {
			for (size_t i = 0; i < obj1.columns(); ++i) {
				obj1(j, i) = obj(i, j);
			}
		}
		return obj1;
	}

	//Функция соеденения матриц*
	const Matrix concatenate(const Matrix obj1, const Matrix obj2) {
		if (obj1.rows()!=obj2.rows())
			throw std::runtime_error("Error, value of rows are not same");

		Matrix obj(obj1.rows(), obj1.columns() + obj2.columns());

		for (size_t i = 0; i < obj.rows(); ++i) {
			for (size_t j = 0; j < obj.columns(); ++j) {
				if (j > obj1.columns() - 1) {
					obj(i, j) = obj2(i, j - obj1.columns());
				}
				else
					obj(i, j) = obj1(i, j);

			}
		}
		return obj;
	}

	//Функция получения строки матрицы*
	const Matrix get_row(const size_t row, const Matrix obj1) {
		if (row >= obj1.rows())
			throw std::runtime_error("Error, invalid row");

		Matrix obj(1, obj1.columns());
		size_t iter = 0;

		for (size_t i = 0; i < obj1.rows() * obj1.columns(); ++i) {
			if (i / obj.columns() == row) {
				obj(iter / obj.columns(), iter % obj.columns()) = obj1(i / obj1.columns(), i % obj1.columns());
				++iter;
			}
		}
		return obj;
	}

	//Функция получения столбца матрицы*
	const Matrix get_column(const size_t column, const Matrix obj1) {
		if (column >= obj1.columns())
			throw std::runtime_error("Error, invalid row");

		Matrix obj(obj1.rows());
		size_t iter = 0;

		for (size_t i = 0; i < obj1.rows() * obj1.columns(); ++i) {
			if (i % obj1.columns() == column) {
				obj(iter / obj.columns(), iter % obj.columns()) = obj1(i / obj1.columns(), i % obj1.columns());
				++iter;
			}
		}
		return obj;
	}

	//Функция получения обратной матрицы*
	const Matrix invert(const Matrix obj) {
		Matrix obj1 = obj, obj_alg_ext = transpose(obj1.mat_alg_ext());

		if (obj1.det() == 0||obj1.columns()!=obj1.rows())
			throw std::runtime_error("Error, cant find invert of this matrix");

		double det_obj = obj1.det();
		obj1 = obj_alg_ext * (1 / det_obj);
		return obj1;
	}

	//Функция возведения матрицы в степень*
	const Matrix power(const Matrix obj, const int val) {
		Matrix obj_old = obj, obj_new = obj;
		if(val>0)
			for (size_t i = 0; i < val - 1; ++i)
				obj_new = obj_new * obj_old;

		if (val == 0) {
			obj_old = invert(obj_old);
			obj_new *= obj_old;
		}

		if(val<0) {
			for (size_t i = 0; i < std::abs(val) - 1; ++i)
				obj_new *= obj_old;
			obj_new = invert(obj_new);
		}
		return obj_new;
	}

	//Функция решения СЛАУ методом Гаусса*
	const Matrix solve(const Matrix obj_mat, const Matrix obj_vct) {
		Matrix obj(obj_mat.columns()), obj_dopol = concatenate(obj_mat, obj_vct), obj_mat_c= obj_mat;
		size_t iter = 0;

		obj_dopol.gauss_forward();

		if(obj_mat_c.rank() < obj_dopol.rank())
			throw std::exception();

		obj_dopol.gauss_backward();

		if (obj_dopol.rank() != obj_mat.columns())
			throw std::exception();

		for (size_t i = 0; i < obj_dopol.columns() * obj_dopol.rows(); ++i) {
			if (i % obj_dopol.columns() == obj_dopol.columns() - 1)
				++i;
			if (obj_dopol(i / obj_dopol.columns(), i % obj_dopol.columns()) != 0 && iter < obj.rows()) {
				obj(iter / obj.columns(), iter % obj.columns()) = \
					obj_dopol(i / obj_dopol.columns(), obj_dopol.columns() - 1)\
					/ obj_dopol(i / obj_dopol.columns(), i % obj_dopol.columns());
				++iter;
			}
		}
		return obj;
	}

	//Определение перегрузки оператора << (Вывод в поток)*
	std::ostream& operator << (std::ostream& out, const Matrix& obj) {
		if (!obj.empty()) {
			std::string str;
			size_t i_row = 0, i_column = 0;
			size_t max_size=0;

			for (size_t i = 0; i < obj.columns() * obj.rows(); ++i) {
				std::ostringstream ss;
				i_row = i / obj.columns();
				i_column = i % obj.columns();
				ss << obj(i_row, i_column);
				str = ss.str();
				if (str.size() > max_size)
					max_size = str.size();
			}

			out << "|";
			for (size_t i = 0; i < obj.columns() * obj.rows(); ++i) {
				if (i % obj.columns() > obj.columns())
					throw std::exception();
				i_row = i / obj.columns();
				i_column = i % obj.columns();	
				if (i_column != obj.columns() - 1) {
					out << std::setw(max_size) \
						<< obj(i_row, i_column) << " ";
				}
				if (i == obj.columns() * obj.rows() - 1) {
					out << std::setw(max_size) \
						<< obj(i_row, i_column) << "|";
				}
				else {
					if (i_column == obj.columns() - 1) {
						out << std::setw(max_size) \
							<< obj(i_row, i_column) << "|\n|";
					}
				}
			}
		}
		return out;
	}
}