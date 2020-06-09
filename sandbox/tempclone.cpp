#include <iostream>

class ab {
    public:
    ab() {
        std::cout<<"Chiamato costruttore"<<std::endl;
    }

    ~ab() {
        std::cout<<"Chiamato distruttore"<<std::endl;
    }

    void ciao() {
        std::cout<<"Ciao"<<std::endl;
    }

    ab* return_self() {
        std::cout<<"ab2"<<this<<std::endl;
        return this;
    }
};

ab& func1() {
    //ab *ab1 = new ab;
    //return *ab1;
    ab &ab1 = *(new ab);
    return ab1;
}

void func2() {

}

int main()
{
    ab &ab2 = func1();
    ab2.ciao();
    ab *ab3 = ab2.return_self();
    ab3->ciao();
    std::cout<<" ab3 "<<ab3<<std::endl;
    delete(ab3);
    ab2.return_self();
    std::cout<<"Dopo tutto"<<std::endl;
    return 0;
    /*
    ab* ab2 = new ab; //func1();
    ab2->ciao();
    ab *ab3 = ab2->return_self();
    ab3->ciao();
    std::cout<<"ab2"<<ab2<<" ab3 "<<ab3<<std::endl;
    delete(ab3);
    std::cout<<"Dopo tutto"<<std::endl;
*/
}