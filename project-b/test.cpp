#include <iostream>

void sleep();

int main()
{
    while(1)
    {
        //system("clear");
        system("date +%r&");
//        sleep(1); 
    }
    return 0;
}

void sleep() {
    int i=0;
    while(i!=1000000){i++;}
}
