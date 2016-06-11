#include<iostream>

#include<fstream>
using namespace std;
int main() 
{
	double a,b;// a is time and b is current
	ofstream out_data("Threshold_variability_data.isf");
	out_data <<"time,0"<<endl;
	for (double a = 100,b = 0; a <=101;a = a + 0.01){
		b = 50 ;
		out_data <<a<<","<<b<<endl;
		}
	for (double a = 110,b = 0; a <=111;a = a + 0.01){
		b = 50 ;
		out_data <<a<<","<<b<<endl;
		}
	for (double a = 120,b = 0; a <=121;a = a + 0.01){
		b = 50 ;
		out_data <<a<<","<<b<<endl;
		}
	for (double a = 130,b = 0; a <=131;a = a + 0.01){
		b = 50 ;
		out_data <<a<<","<<b<<endl;
		}
	return 0;
}

