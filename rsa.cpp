#include <string>
#include <iostream>
#include<cmath>
#define BYTESIZE 256/8
void print(unsigned char* a){
    for (int i = 0; i < BYTESIZE;i++){
        std::cout << a[i]-0 << " ";
    }
    std::cout << std::endl;
}

bool char_self_add(unsigned char& a, unsigned char b, bool add){
    bool overflow = false;
    unsigned char tmp = 0;
    while(b){
        tmp = b;
        b &= a;
        if(b >= 0x80){
            overflow = true;
        }
        b <<= 1;
        a ^= tmp;
    }
    if(add){
        b = 1;
        while(b){
            tmp = b;
            b &= a;
            if(b >= 0x80){
                overflow = true;
            }
            b <<=  1;
            a ^= tmp;
        }
    }
    return overflow;
}

bool add_bit(unsigned char* a,unsigned char *b){
    bool add = false;
    for (int i = 0; i < BYTESIZE; i++){
        add = char_self_add(a[i], b[i], add);
    }
    return add;
}
bool minus_bit(unsigned char* a,unsigned char *b){
    bool add = true;
    for (int i = 0; i < BYTESIZE; i++){
        add = char_self_add(a[i], ~b[i], add);
    }
    return add;
}

bool compare(const unsigned char* data, const unsigned char* bytes){
    for (int i = BYTESIZE-1; i >= 0; i--){
        if(data[i]>bytes[i]){
            return true;
        }else if(data[i]<bytes[i]){
            return false;
        }
    }
    return false;
}
bool equal(const unsigned char* data, const unsigned char* bytes){
    for (int i = 0; i < BYTESIZE; i++){
        if(data[i]!=bytes[i]){
            return false;
        }
    }
    return true;
}
bool int2bi(int a, unsigned char* result){
    int i = 0;
    while (i<BYTESIZE*8 &&  a)
    {
        int first = i / 8;
        int second = i % 8;
        if(a%2 == 1){
            result[first] |= (1 << second);
        }
        i += 1;
        a /= 2;

    }

    return true;
}

int str2char(char* result, unsigned char* text){
    int first = 0;
    int second = 0;
    int size = 0;
    int value = 0;
    while (result[first])
    {
        if(size==8/4){
            text[second] = value;
            size = 0;
            second += 1;
            value = 0;
        }
        if(result[first]<='9'&& result[first]>='0'){
            value += (pow(16, size)*(result[first]-'0'));
        }else if(result[first]<='f'&& result[first]>='a'){
            value += (pow(16, size)*(result[first]-'a'+10));
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

int char2str(unsigned char *encode_text, int size, char* result){
    int index = 0;
    int answer = 0;
    for (int i = 0; i < size; i++){
        index = encode_text[i];
        int j = 0;
        while (index)
        {
            answer = index % 16;
            if (answer <= 9)
            {
                result[i * 8 / 4 + j] = '0'+answer;
            }
            else
            {
                result[i * 8 / 4 + j] = 'a'+(answer-10);
            }
            index = index / 16;
            j += 1;
        }
        for (; j < 8 / 4;j++){
            result[i * 8 / 4 + j] = 'g';
        }
    }
    result[8 / 4 * size + 1] = 0;
    return 8 / 4 * size + 1;
}
class BytesBase{
private:
    unsigned char _base[2*BYTESIZE*BYTESIZE*8] = {0};
    unsigned char _data[BYTESIZE];
    int _large_index = BYTESIZE*8+1;

public:
    BytesBase(){}
    BytesBase(const unsigned char *bytes)
    {
        for (int i = 0; i<BYTESIZE;i++){
            _data[i] = bytes[i];
        }
        int value= 1;
        unsigned char tmp_base[BYTESIZE] = {0};
        unsigned char bias[BYTESIZE] = {0};
        bias[0] = 1;
        tmp_base[0] = 1;
        for (int i = 0; i < BYTESIZE * 2; i++){
            for (int j = 0; j < 8; j++){

                for (int t = 0; t < BYTESIZE; t++){
                    _base[(i * 8 + j) * BYTESIZE + t] = tmp_base[t];
                }
                if(i<BYTESIZE){
                    if((_base[(i * 8 + j) * BYTESIZE+i]&(1<<j))==0){
                        if(_large_index>i * 8 + j){
                            _large_index = i * 8 + j;
                        }
                    }
                }

                value = (value * 2) % 2123;
                bool add = add_bit(tmp_base, tmp_base);
                while (add){
                    add = add_bit(tmp_base, _base + (BYTESIZE * 8) * BYTESIZE);
                }
                if(i < BYTESIZE){
                    int h = i;
                    int k = j + 1;
                    if (k == 8){
                        k = 0;
                        h += 1;
                    }
                    for (int s = 0; s < BYTESIZE; s++){
                        if(s==h){
                            tmp_base[s]=((1<<k) - 1);
                        }else if(s<h) {
                            tmp_base[s]=((1<<8) - 1);
                        }else{
                            tmp_base[s] = 0;
                        }
                    }
                    get_residual(tmp_base, tmp_base);
                    add_bit(tmp_base, bias);
                }else{
                    get_residual(tmp_base, tmp_base);
                }
            }
        }
    }

    void    get_residual(const unsigned char* bytes, unsigned char* ans){
        if(compare(_data, bytes)){
            for (int i = 0; i < BYTESIZE; i++){
                ans[i] = bytes[i];
            }
            return;
        }
        unsigned char tmp_ans[BYTESIZE] = {0};
        bool change =false;
        for (int i = 0; i < BYTESIZE; i++){
            for (int j = 0; j < 8;j++){
                if(bytes[i] & (1 << j)){
                    add_bit(tmp_ans, _base+((i*8+j)*BYTESIZE));

                    if(8 * i + j >= _large_index){
                        change = true;
                    }
                }
            }
        }
        if (!change){
            minus_bit(tmp_ans, _data);
        }
        if(equal(_data, tmp_ans)){
            for (int i = 0; i < BYTESIZE; i++){
                ans[i] = 0;
            }
            return;
        }else if(compare(_data, tmp_ans)){
            for (int i = 0; i < BYTESIZE; i++){
                ans[i] = tmp_ans[i];
            }
            return;
        }else{
            get_residual(tmp_ans, ans);
            return;
        }
    }

    void exp(const unsigned char* exp_base, unsigned char* exp, unsigned char*  ans){
        int index = 0;
        unsigned char tmp_ans[BYTESIZE];
        for (int i = 0; i<BYTESIZE;i++){
            tmp_ans[i] = exp_base[i];
        }
        for (int i = 1; i<BYTESIZE;i++){
            ans[i] =0;
        }
        int first = 0;
        int second = 0;
        bool add = false;
        ans[0] = 1;
        while (index < BYTESIZE)
        {
            first = index / 8;
            second = index % 8;
            if(exp[first]&1<<second){
                multiply(ans, tmp_ans, ans);
            }
            index += 1;
            multiply(tmp_ans, tmp_ans, tmp_ans);
        }
        get_residual(ans, ans);
    }

    void multiply(unsigned char* a,  unsigned char* b,  unsigned char* ans){
        unsigned char tmp_ans[BYTESIZE] = {0};
        for (int i = 0; i < BYTESIZE;i++){
            for (int j = 0; j < 8;j++){
                if(a[i]&(1<<j)){
                    for (int s = 0; s < BYTESIZE; s++){
                        for (int t = 0; t < 8;t++){
                            if(b[s]&(1<<t)){
                                bool add = add_bit(tmp_ans, _base + (i * 8 + j + s * 8 + t) * BYTESIZE);
                                while(add){
                                    add = add_bit(tmp_ans, _base + 8*BYTESIZE * BYTESIZE);
                                }
                            }
                        }
                    }
                }
            }
        }
        get_residual(tmp_ans, ans);
    }
};

class Encrpyt{
    private:
        BytesBase _base_n;
        unsigned char _value_e[BYTESIZE]={0};
        unsigned char _key_d[BYTESIZE]={0};
        int _size;

    public:
        Encrpyt(int p, int q, int e){
            unsigned char tmp[BYTESIZE] = {0};
            int2bi(p * q, tmp);
            _base_n = BytesBase(tmp);
            for (int i = 0;i<BYTESIZE; i++){
                if((tmp[i]&((8<<1)-1))){
                    _size = i;
                }
            }
                int2bi(e, _value_e);


            unsigned char tmp1[BYTESIZE] = {0};
            int2bi((p-1) * (q-1), tmp1);
            BytesBase euler_base = BytesBase(tmp1);

            unsigned char bias[BYTESIZE] = {0};
            unsigned char tmp2[BYTESIZE] = {0};
            bias[0] = 1;
            bool add = true;
            while (compare(tmp1, _key_d)&& add)
            {
                euler_base.get_residual(tmp2, tmp2);
                if (equal(tmp2, bias)){
                    add = false;
                }else{
                    add_bit(_key_d, bias);
                    add_bit(tmp2, _value_e);
                }
                
            }
        }
        void encrp_data(unsigned char* text, unsigned char* result){
            _base_n.exp(text, _value_e, result);
        }
        void decrp_data(unsigned char* result, unsigned char* text){
            _base_n.exp(result, _key_d, text);
        }

        int encrp(const char* a, unsigned char* encode_text){
            unsigned char text[BYTESIZE]={0};
            unsigned char seg_result[BYTESIZE]={0};
            int encode_size = (sizeof(a)-1)/_size+1;
            int text_index = 0;
            int encode_index = 0;
            int i=0;
            while(a[i]){
                if(text_index==_size){
                    encrp_data(text, seg_result);
                    for (int i = 0; i <_size+1;i++){
                        encode_text[encode_index] = seg_result[i];
                        encode_index += 1;
                    }
                    text_index = 0;
                }
                text[text_index] = a[i];
                text_index += 1;
                i+=1;
            }
            if(text_index==_size){
                encrp_data(text, seg_result);
                for (int i = 0; i < _size+1;i++){
                    encode_text[encode_index] = seg_result[i];
                    encode_index += 1;
                }
            }else{
                for (int i = text_index; i<_size;i++){
                    text[text_index] = 0;
                }
                encrp_data(text, seg_result);
                for (int i = 0; i < get_size()+1;i++){
                    encode_text[encode_index] = seg_result[i];
                    encode_index += 1;
                }
            }
            return encode_index;

        }

        void decrp(unsigned char *encode_text, int size,  unsigned char* decode_text)
        {
            int decode_index = 0;
            int text_index = 0;
            unsigned char seg_encode[BYTESIZE]={0};
            unsigned char seg_result[BYTESIZE]={0};
            for (int i = 0; i < size; i++){
                if(text_index == _size +1){
                    decrp_data(seg_encode, seg_result);
                    for (int j = 0; j < _size; j++){
                       
                        decode_text[decode_index] = seg_result[j];
                        decode_index += 1;
                    }
                    text_index = 0;
                }
                seg_encode[text_index] = encode_text[i];
                text_index += 1;
            }
            if(text_index == _size +1){
                decrp_data(seg_encode, seg_result);
                for (int j = 0; j < get_size();j++){
                    decode_text[decode_index] = seg_result[j];
                    decode_index += 1;
                }
                text_index = 0;
            }
        }
        void decrp(char* o,  unsigned char *decode_text1 ){
            unsigned char* encode_test1=new unsigned char[sizeof(o)/(8/4)+1];
            int size = str2char(o, encode_test1);
            decrp(encode_test1, size , decode_text1);
            //delete[] encode_test1;
        }
         void encrp(const char* o,  char* result ){
            unsigned char* encode_text = new unsigned char[((sizeof(o)-1)/_size+1)*(_size+1)];
            int size=encrp(o, encode_text);
            size = char2str(encode_text, size, result);
            //delete[] encode_text;
        }
        int get_size(){
            return _size;
        }
};




int main(){

    std::string a = "你好,你叫什么名字";
    Encrpyt enc(647,653,659);
    int encode_size = ((a.size() - 1) / enc.get_size() + 1)*(enc.get_size()+1)*2+1;
    char result[encode_size];
    enc.encrp(a.c_str(), result );
    std::cout<<result<<std::endl;
    
    unsigned char decode_text[(sizeof(result)-1)*enc.get_size()/2/(enc.get_size()+1)+1];
    enc.decrp(result, decode_text);
    std::cout << (char*)decode_text << std::endl;
    return 0;
}
