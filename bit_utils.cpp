
#include "bit_base.cpp"
template<typename T>
bool int2bi(int a, T* result){
    int i = 0;
    while (i<BYTESIZE*8*sizeof(T) &&  a)
    {
        int first = (i/(sizeof(T)*8));
        int second = (i&(8*sizeof(T)-1));
        if((a&1) == 1){
            result[first] |= (((T)1) << second);
        }
        i += 1;
        a /= 2;

    }

    return true;
}

template<typename T>
int str2char(char* result,int char_size, T* text){
    int first = 0;
    int second = 0;
    int size = 0;
    T value = 0;
    for (int i = 0; i < char_size;i++)
    {
        if(size==8*sizeof(T)/4){
            text[second] = value;
            size = 0;
            second += 1;
            value = 0;
        }
        if(result[first]<='9'&& result[first]>='0'){
            value += (pow(((T)16), size)*((T)(result[first]-'0')));
        }else if(result[first]<='f'&& result[first]>='a'){
            value += (pow(((T)16), size)*((T)(result[first]-'a'+10)));
        }else if(result[first]!='g'){
            return 0;
        }
        size += 1;
        first += 1;
    }
    text[second] = value;
    size = 0;
    second += 1;
    value = 0;
    return second;
}

template<typename T>
int char2str(T *encode_text, int size, char* result){
    T index = 0;
    int answer = 0;
    for (T i = 0; i < size; i++){
        index = encode_text[i];
        int j = 0;
        while (index)
        {
            answer = (index & 15);
            if (answer <= 9)
            {
                result[i * 8*sizeof(T) / 4 + j] = '0'+answer;
            }
            else
            {
                result[i * 8*sizeof(T) / 4 + j] = 'a'+(answer-10);
            }
            index >>=4;
            j += 1;
        }
        for (; j < 8*sizeof(T) / 4;j++){
            result[i * 8*sizeof(T) / 4 + j] = 'g';
        }
    }
    result[8 / 4*sizeof(T) * size] = 0;
    return 8 / 4*sizeof(T) * size + 1;
}

