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
bool equal_zero(const unsigned char* data){
    for (int i = 0; i < BYTESIZE; i++){
        if(data[i]!=0){
            return false;
        }
    }
    return true;
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
int log2(unsigned char* a){
     int size = 0;
    int start = BYTESIZE*8-1;
    int first =0;
    int second =0;
    while (size==0&&start>0){
        first=start&7;
        second = start>>3;
        if(a[second]&(1<<first)){
            return start;
        }
        start-=1;
    }
    return -1;
}
void left_bit(unsigned char* a,int size){
    int pre_fist=(BYTESIZE-1-(size>>3));
    int pre_bias =(size&7);
    unsigned char pre=0;
    for(int i=BYTESIZE-1;i>=0;i--){
        if(pre_fist<0){
            a[i]=0;
        }else if(pre_fist==0){
            a[i]=(a[pre_fist]<<pre_bias);
        }else{
            a[i]=(a[pre_fist]<<pre_bias)+(a[pre_fist-1]>>(8-pre_bias));
        }
        pre_fist-=1;
    }
}
void right_bit(unsigned char* a,int size){
    int pre_fist=(size>>3);
    int pre_bias =(size&7);
    unsigned char pre=0;
    for(int i=0;i<BYTESIZE;i++){
        if(pre_fist>BYTESIZE){
            a[i]=0;
        }else if(pre_fist==BYTESIZE-1){
            a[i]=(a[pre_fist]>>(pre_bias));
        }else{
            a[i]=(a[pre_fist]>>(pre_bias))+(a[pre_fist+1]<<(8-pre_bias));
        }
        pre_fist+=1;
    }
}

bool divide_bit(unsigned char* a,unsigned char *b,unsigned char *redis){
    unsigned char ans[BYTESIZE] ={0};
    unsigned char bias[BYTESIZE];
    for(int i=0;i<BYTESIZE;i++){
        bias[i]=b[i];
    }
    int size = log2(a)-log2(b);
    left_bit(bias,size);

    while(!compare(b,a)){
        
        if(!compare(bias,a)){
            ans[size>>3]+=(1<<(size&7));
            minus_bit(a,bias);
        }
        right_bit(bias,1);
        size-=1;
    }
    for(int i=0;i<BYTESIZE;i++){
       redis[i]=a[i];
       a[i]=ans[i];
    }
    
}

bool 
multiply_bit(unsigned char* a,unsigned char *b, unsigned char*ans){
    bool a_flat[BYTESIZE*8];
    bool b_flat[BYTESIZE*8];
    for(int i=0;i<BYTESIZE*8;i++){
        int first =i>>3;
        int second =(1<<(i&7));
        if(a[first]&second){
            a_flat[i]=true;
        }else{
            a_flat[i]=false;
        }
        if(b[first]&second){
            b_flat[i]=true;
        }else{
            b_flat[i]=false;
        }

    }
    int first=0;
    int second=0;
    int number = 1;
    unsigned char bias[BYTESIZE]={0};
    while(first<BYTESIZE){
        bias[first]=(1<<second);
        for(int i=0;i<number;i++){
            if(a_flat[i]&&b_flat[number-1-i]){
                if(add_bit(ans, bias)){
                    return true;
                }
            }
        }
        second+=1;
        number+=1;
        if(second==8){
            bias[first]=0;
            second = 0;
            first+=1;
            if(first<BYTESIZE){
                bias[first]=1;
            }
        }
    }
    int start =1;
    while (number<2*BYTESIZE*8) {
        for(int i=start;i<number;i++){
            if(a_flat[i]&&b_flat[number-1-i]){
                return true;
            }
        }
        number+=1;
    }
    
    return false;
}
bool ex_gcd(unsigned char* a, unsigned char* b, unsigned char* x, unsigned char* y, bool& x_flag, bool& y_flag){
    if(equal_zero(b)){
        for(int i=0;i<BYTESIZE;i++){
            if(i==0){
                if(a[i]!=1){
                    return false; 
                }
                x[i]=1;
                y[i]=0;
            }else{
                if(a[i]!=0){
                    return false; 
                }
                x[i]=0;
                y[i]=0;
            }
        }
        x_flag=true;
        y_flag =true;
        return true;
    }else{
        unsigned char next_a[BYTESIZE];
        unsigned char next_b[BYTESIZE];
        unsigned char resi[BYTESIZE];
        for(int i=0;i<BYTESIZE;i++){
            next_a[i]=a[i];
            next_b[i]=b[i];
        }
        divide_bit(next_a,next_b,resi);        
        if(!ex_gcd(next_b,resi,y,x,y_flag,x_flag)){
            return false;
        }
        for(int i=0;i<BYTESIZE;i++){
            next_b[i]=0;
        }
        multiply_bit(next_a,x,next_b);
        if(y_flag&&x_flag){
            if(!compare(next_b,y)){
                minus_bit(y,next_b);
            }else{
                minus_bit(next_b,y);
                for(int i=0;i<BYTESIZE;i++){
                    y[i]=next_b[i];
                }
                y_flag=false;
            }
        }else if(y_flag&&!x_flag){
             if(add_bit(y,next_b)){
                 return false;
             }
        }else if(x_flag&&!y_flag){
             if(add_bit(y,next_b)){
                 return false;
             }
        }else{
            if(!compare(next_b,y)){
                minus_bit(y,next_b);
            }else{
                minus_bit(next_b,y);
                for(int i=0;i<BYTESIZE;i++){
                    y[i]=next_b[i];
                }
                y_flag=true;
            }
        }
        return true;
    }
}

bool inv(unsigned char* t, unsigned char* p, unsigned char* ans){
    unsigned char x[BYTESIZE];
    bool x_flag;
    bool y_flag;
    if(!ex_gcd(t,p,ans,x,x_flag,y_flag)){
        return false;
    };
    if(!x_flag){
        divide_bit(ans,p,x);
        for(int i=0;i<BYTESIZE;i++){
            ans[i]=p[i];
        }
        minus_bit(ans,x);
    }

    return true;
}


bool int2bi(int a, unsigned char* result){
    int i = 0;
    while (i<BYTESIZE*8 &&  a)
    {
        int first = (i>>3);
        int second = (i&7);
        if((a&1) == 1){
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
            answer = (index & 15);
            if (answer <= 9)
            {
                result[i * 8 / 4 + j] = '0'+answer;
            }
            else
            {
                result[i * 8 / 4 + j] = 'a'+(answer-10);
            }
            index >>=4;
            j += 1;
        }
        for (; j < 8 / 4;j++){
            result[i * 8 / 4 + j] = 'g';
        }
    }
    result[8 / 4 * size] = 0;
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
    void print_bit(){
        for(int i=0;i<BYTESIZE*2*8;i++){
            print(_base+i*BYTESIZE);
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
        int size = log2(exp)+1;
        while (index < size)
        {
            first = (index>>3);
            second = (index&7);
            if(exp[first]&1<<second){
                multiply(ans, tmp_ans, ans);
            }
            index += 1;
            multiply(tmp_ans, tmp_ans, tmp_ans);
        }
        get_residual(ans, ans);
    }
    
    void  multiply(unsigned char* a,unsigned char *b, unsigned char*ans){
        bool a_flat[BYTESIZE*8];
        bool b_flat[BYTESIZE*8];
        for(int i=0;i<BYTESIZE*8;i++){
            int first =(i>>3);
            int second =(1<<(i&7));
            if(a[first]&second){
                a_flat[i]=true;
            }else{
                a_flat[i]=false;
            }
            if(b[first]&second){
                b_flat[i]=true;
            }else{
                b_flat[i]=false;
            }

        }
        int number = 0;
        unsigned char tmp_ans[BYTESIZE] = {0};
        while(number<2*BYTESIZE*8){
            for(int i=0;i<number+1;i++){
                if(i<BYTESIZE*8&&number-i<BYTESIZE*8){
                    if(a_flat[i]&&b_flat[number-i]){
                        bool add = add_bit(tmp_ans, _base+number*BYTESIZE);
                        while(add){
                            add = add_bit(tmp_ans, _base + 8*BYTESIZE * BYTESIZE);
                        }
                    }
                }
                
            }
            number+=1;

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
        bool _ok;

    public:
        Encrpyt(int p, int q, int e){
            unsigned char tmp[BYTESIZE] = {0};
            unsigned char tmp2[BYTESIZE] = {0};
            unsigned char tmp3[BYTESIZE] = {0};
            int2bi(p,tmp2);
            int2bi(q,tmp3);
            multiply_bit(tmp2,tmp3,tmp);
            _base_n = BytesBase(tmp);
            for (int i = 0;i<BYTESIZE; i++){
                if((tmp[i]&((8<<1)-1))){
                    _size = i;
                }
            }
            unsigned char bias[BYTESIZE]={0};
            bias[0]=1;
            minus_bit(tmp2,bias);
            minus_bit(tmp3,bias);
            int2bi(e, _value_e);
            unsigned char tmp1[BYTESIZE] = {0};
             multiply_bit(tmp2,tmp3,tmp1);
            _ok=inv(_value_e,tmp1,_key_d);
            if(!_ok){
                std::cout<<"init failed"<<std::endl;
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
    Encrpyt enc(500881,500887,500891);
    int encode_size = ((a.size() - 1) / enc.get_size() + 1)*(enc.get_size()+1)*2+1;
    char result[encode_size];
    enc.encrp(a.c_str(), result );
    std::cout<<result<<std::endl;
    
    unsigned char decode_text[(sizeof(result)-1)*enc.get_size()/2/(enc.get_size()+1)+1];
    enc.decrp(result, decode_text);
    std::cout << (char*)decode_text << std::endl;
    return 0;
}
