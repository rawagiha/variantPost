#include <iostream>
#include <string>

#include "swlib.h"

int main() {
    
    std::string read1 =  "GTTAAAAAAGTATTTTCAGGGCTTGGAGCGGTGG";
    std::string qread1 = "AAAFFJJJJJJJJJJJJAF<FJAJJJJJJJJJJB";
    std::string read2 = "AAAAGTATTGTCAGAGCTTGGAGCGGTGGCTCTTG";
    //std::string read2 = "GAGCGGTGGAGTATTGTCTCTTGAAACAGAGCTTG ";
    std::string qread2 = "J-F---7FJFA-AAJF-<-7<<-AJF----ADFGA";
     
    //std::string read3 =  "GCTTGGAGCGGTGGCTCTTGCCTGTAATGCCAGCACTTTGGGAGGCCAAGGCAGGCAGA";
    
    std::string read3 =    "GCTTGGNGCGGTGGCTCTTGCCTGTAATGCCAGCACTTTGGGAGGCCAAGGCAGGCAGA";
    std::string qread3 = "-AAJF-<-7<<-AJF----AAFFJJJJJJJJJJAAAFFJJJJAAJF-<-7<<-AJF---";  
    
    std::string read4 =  "GATGACTAGTATTA";
    std::string qread4 = "AAAAAAAAAAAAAA";

    std::cout << read1.size() << ", " << qread1.size() << std::endl;
    std::cout << read2.size() << ", " << qread2.size() << std::endl;
    std::cout << read3.size() << ", " << qread3.size() << std::endl;

    std::vector<std::string> v1 = {read1, qread1};
    std::vector<std::string> v2 = {read2, qread2};
    std::vector<std::string> v3 = {read3, qread3};
    std::vector<std::string> v4 = {read4, qread4};

    std::vector<std::vector<std::string>> reads = {v1, v4};
    std::string s = sw::flatten_reads(v3, reads);
   
    double d = 10004.3;
    
    double dd = 81.9;
    double ddd = 40.3;
    std::cout << ++dd  << std::endl;
    std::cout << static_cast<char>(d) << std::endl;
    std::cout << static_cast<char>(dd) << std::endl;
    std::cout << static_cast<char>(ddd) << std::endl;

    
}
