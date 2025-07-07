#pragma once
#include <iostream>	
#include <string>
#include <sstream>
#include <iomanip>
namespace linalg
{
	class Matrix {
	public://Конструкторы и деструктор

		Matrix() = default;
		Matrix(int rows);
		Matrix(int rows, int columns);
		Matrix(const Matrix& obj);
		Matrix(Matrix&& obj) noexcept;
		Matrix(std::initializer_list<double> lst);
		Matrix(std::initializer_list<std::initializer_list<double>> lst);
		~Matrix() { delete[] m_ptr; }

	public://Методы-функции

		const size_t rows() const noexcept { return m_rows; }
		const size_t columns() const noexcept{ return m_columns; }
		const bool empty() const noexcept;
		void reshape(const size_t rows, const size_t columns);
		const double norm() const;
		const double trace() const;
		const double det() const;
		const size_t rank() const;
		Matrix& gauss_forward();
		Matrix& gauss_backward();

		//Дополнительные методы
		Matrix& del_row(const size_t row);
		Matrix& del_column(const size_t column);
		const double minor_ext(const size_t row, const size_t column);
		Matrix mat_minor_ext();
		Matrix mat_alg_ext();
		Matrix& swap_rows(const size_t row1, const size_t row2);
		const double sum_el();

	public://Операторы

		Matrix& operator = (const Matrix& obj);
		Matrix& operator = (Matrix&& obj) noexcept ;
		const double operator () (const size_t row, const size_t column) const;
		double& operator () (const size_t row, const size_t column);
		Matrix& operator += (const Matrix& obj);
		Matrix& operator *= (const Matrix& obj);
		Matrix& operator -= (const Matrix& obj);
		const Matrix operator + (const Matrix& obj) const ;
		const Matrix operator * (const Matrix& obj) const ;
		const Matrix operator - (const Matrix& obj) const ;
		Matrix& operator *= (const double val);
		const Matrix operator * (const double val) const ;

	private://Ресурсы класса

		size_t m_rows=0, m_columns=0;
		double* m_ptr=nullptr;
	};

	//Функции требующие глобального объявления 
	const Matrix transpose(const Matrix obj);
	const Matrix concatenate(const Matrix obj1, const Matrix obj2);
	const Matrix get_row(const size_t row, const Matrix obj1);
	const Matrix get_column(const size_t column, const Matrix obj1);
	const Matrix invert(const Matrix obj);
	const Matrix power(const Matrix obj, const int val);
	const Matrix solve(const Matrix obj_mat, const Matrix obj_vct);
	//Операторы требующие глобального объявления
	std::ostream& operator << (std::ostream& out, const Matrix& obj);
	const bool operator == (const linalg::Matrix& obj1, const linalg::Matrix& obj2);
	const bool operator != (const linalg::Matrix& obj1, const linalg::Matrix& obj2);
	const linalg::Matrix operator * (const double val, const linalg::Matrix& obj);
}
