#include <iostream>

class SizeException {
    const char* const data;
public:
    SizeException(const char* const msg = 0) : data(msg) {}
};

void Err() {
    throw SizeException("La mère de Nicolas est une p*te")
}

int main() {
    Err();
}