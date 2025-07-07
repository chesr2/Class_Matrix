#include "Test.h"
namespace tst {
	void test1(linalg::Matrix& obj) {
		size_t m, n;
		std::cout << "Your matrix is - \n" << obj << "\n\n";
		std::cout << "Value of rows: " << obj.rows()\
			<< "\n" << "Value of columns: " << obj.columns() << "\n\n";
		std::cout << "Is your matrix empty or not? Answer - " << obj.empty() << "\n\n";
		std::cout << "Reshape your Matrix to (Enter value of rows and columns):";
		std::cin >> m;
		std::cout << "x";
		std::cin >> n;
		obj.reshape(m, n);
		std::cout << "Result:\n" << obj << "\nEND OF TEST 1\n\n";
	}

	void test2() {
		linalg::Matrix obj1;
		std::cout << "linalg::Matrix obj1; - \n" << obj1 << "\n\n";
		linalg::Matrix obj2(4);
		std::cout << "linalg::Matrix obj2(4); - \n" << obj2 << "\n\n";
		linalg::Matrix obj3(3, 2);
		std::cout << "linalg::Matrix obj3(3,2); - \n" << obj3 << "\n\n";
		linalg::Matrix obj4(obj3);
		std::cout << "linalg::Matrix obj4(obj3); - \n" << obj4 << "\n\n";
		linalg::Matrix obj5(std::move(obj4));
		std::cout << "linalg::Matrix obj5(std::move(obj4)); - \n" << "linalg::Matrix obj5(obj4); - \n" << "obj5 - \n" << obj5 << "\nobj4 - \n" << obj4 << "\n\n";
		linalg::Matrix obj6 = { {1,2,3},{4,5,6},{7,8,9} };
		std::cout << "linalg::Matrix obj6 = { {1,2,3},{4,5,6},{7,8,9} }; - \n" << obj6 << "\n\n";
		linalg::Matrix obj7 = { {1,2,3,4,5,6} };
		std::cout << "linalg::Matrix obj7 = { {1,2,3,4,5,6} }; - \n" << obj7 << "\n\n";
		linalg::Matrix obj8 = { 1,2,3,4,5,6 };
		std::cout << "linalg::Matrix obj8 = { 1,2,3,4 }; - \n" << obj8 << "\n\n";
		linalg::Matrix obj9 = { {1},{2},{3},{4},{5},{6} };
		std::cout << "linalg::Matrix obj9= { {1},{2},{3},{4} }; - \n" << obj9 << "\nEND OF TEST 2\n";
	}

	void test3(linalg::Matrix& obj1, linalg::Matrix& obj2) {
		obj1 = obj2;
		std::cout << "obj1 = obj2;(Copy) -\nobj1- \n" << obj1 << "\nobj2 - \n" << obj2 << "\n";
		obj1 = linalg::Matrix{ 1, 2, 3, 4, 5, 6 };
		std::cout << "obj1 = linalg::Matrix{ 1, 2, 3, 4, 5, 6 };(Replace) -\nobj1- \n" << obj1 << "\n";
		std::cout << "What value have obj2(0, 2)? - " << obj2(0, 2) << "\n";
		obj2(0, 2) = 7;
		std::cout << "obj2(0, 2) = 7 - \n" << obj2 << "\n";
		const linalg::Matrix m = { {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0} };
		double val = m(0, 2);
		//m(0, 2) = 7.0;
		std::cout << "const linalg::Matrix m = { {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0} }; - m(0, 2) = 7 - MISTAKE \n" << "\nEND OF TEST 3\n";
	}

	void test4(linalg::Matrix& obj) {
		linalg::Matrix obj1 = { {10,0.1,1},{2.5,34,5},{12,3,0.8} }, obj2 = { {132454,1},{1,0.32848} },\
			obj3 = { 123456,12345,1234,123,12,1 };
		std::cout << "Your matrix - \n" << obj << "\n";
		std::cout << "What can do << else? -\n" << obj1 << "\n\n" << obj2\
			<< "\n\n" << obj3 << "\nEND OF TEST 4";
	}

	void test6() {

	}

	void test7() {

	}

	void test8() {

	}
}
