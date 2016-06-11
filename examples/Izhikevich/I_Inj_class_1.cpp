#include<iostream>

#include<fstream>
using namespace std;
int main() 
{
	double a,b;// a is time and b is current
	ofstream out_data("Class_1_excitable_data.isf");
	out_data <<"time,0"<<endl;
	for (double a = 0,b = 0; a <=1000;a = a + 0.01){
		b = 0.05*a ;
		out_data <<a<<","<<b<<endl;
		}
	return 0;
}

