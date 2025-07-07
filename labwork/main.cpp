#include "matrix.h"
#include "Test.h"
int main() {
	try {
		linalg::Matrix obj1 = { {12,3,4},{2,334,1},{6,3,2} }, obj2 = { {1,2,3} }, obj(obj1);
		std::cout << obj1;
	}
	catch (std::runtime_error e) {
		std::cerr << e.what();
	}
	catch(std::exception()) {
		std::cerr << "Some Exception";
	}
	catch (...){
		std::cout << "something bad";
	}
}