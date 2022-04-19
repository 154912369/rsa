#include <string>
#include <iostream>
#include<cmath>
#include "bit_base.cpp"
#include "bit_utils.cpp"

template<class T>
class BytesBase{
private:
    T _base[2*BYTESIZE*BYTESIZE*8*sizeof(T)] = {0};
    T _data[BYTESIZE];
    int _large_index = BYTESIZE*8*sizeof(T)+1;

public:
    BytesBase(){}
    BytesBase(const T *bytes)
    {
        for (int i = 0; i<BYTESIZE;i++){
            _data[i] = bytes[i];
        }
        T tmp_base[BYTESIZE] = {0};
        T bias[BYTESIZE] = {0};
        bias[0] = 1;
        tmp_base[0] = 1;
        for (int i = 0; i < BYTESIZE * 2; i++){
            for (int j = 0; j < 8*sizeof(T); j++){

                for (int t = 0; t < BYTESIZE; t++){
                    _base[(i * 8*sizeof(T) + j) * BYTESIZE + t] = tmp_base[t];
                }
                if(i<BYTESIZE){
                    if(((_base[(i * 8*sizeof(T) + j) * BYTESIZE+i]>>j)&1)==0){
                        if(_large_index>i * 8*sizeof(T) + j){
                            _large_index = i * 8*sizeof(T) + j;
                        }
                    }
                }

                bool add = add_bit(tmp_base, tmp_base);
                while (add){
                    add = add_bit(tmp_base, _base + (BYTESIZE * 8*sizeof(T)) * BYTESIZE);
                }
                if(i < BYTESIZE){
                    int h = i;
                    int k = j + 1;
                    if (k == 8*sizeof(T)){
                        k = 0;
                        h += 1;
                    }
                    for (int s = 0; s < BYTESIZE; s++){
                        if(s==h){
                            tmp_base[s]=((((T)1)<<k) - 1);
                        }else if(s<h) {
                            tmp_base[s]=(~((T)0));
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
        for(int i=0;i<BYTESIZE*2*8*sizeof(T);i++){
            print(_base+i*BYTESIZE);
        }
    }


    void    get_residual(const T* bytes, T* ans){
        if(compare(_data, bytes)){
            for (int i = 0; i < BYTESIZE; i++){
                ans[i] = bytes[i];
            }
            return;
        }
        T tmp_ans[BYTESIZE] = {0};
        bool change =false;
        for (int i = 0; i < BYTESIZE; i++){
            for (int j = 0; j < 8*sizeof(T);j++){
                if(bytes[i] & (((T)1) << j)){
                    add_bit(tmp_ans, _base+((i*8*sizeof(T)+j)*BYTESIZE));

                    if(8 * i*sizeof(T) + j >= _large_index){
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

    void exp(const T* exp_base, T* exp, T*  ans){
        int index = 0;
        T tmp_ans[BYTESIZE];
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
            first = (index/(sizeof(T)*8));
            second = (index&(8*sizeof(T)-1));
            if(exp[first]&(((T)1)<<second)){
                multiply(ans, tmp_ans, ans);
            }
            index += 1;
            multiply(tmp_ans, tmp_ans, tmp_ans);
        }
        get_residual(ans, ans);
    }
    
    void  multiply(T* a, T *b, T* ans){
        bool a_flat[BYTESIZE*8*sizeof(T)];
        bool b_flat[BYTESIZE*8*sizeof(T)];
        for(int i=0;i<BYTESIZE*8*sizeof(T);i++){
            int first =(i/(sizeof(T)*8));
            T second =(((T)1)<<(i&(8*sizeof(T)-1)));
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
        T tmp_ans[BYTESIZE] = {0};
        while(number<2*BYTESIZE*8*sizeof(T)){
            for(int i=0;i<number+1;i++){
                if(i<BYTESIZE*8*sizeof(T)&&number-i<BYTESIZE*8*sizeof(T)){
                    if(a_flat[i]&&b_flat[number-i]){
                        bool add = add_bit(tmp_ans, _base+number*BYTESIZE);
                        while(add){
                            add = add_bit(tmp_ans, _base + 8*BYTESIZE * BYTESIZE*sizeof(T));
                        }
                    }
                }
                
            }
            number+=1;

        }
        get_residual(tmp_ans, ans);
    }
};

template<typename T>
class Encrpyt{
    private:
        BytesBase<T> _base_n;
        T _value_e[BYTESIZE]={0};
        T _key_d[BYTESIZE]={0};
        int _size;
        bool _ok;

    public:
        Encrpyt(int p, int q, int e){
            T tmp[BYTESIZE] = {0};
            T tmp2[BYTESIZE] = {0};
            T tmp3[BYTESIZE] = {0};
            int2bi(p,tmp2);
            int2bi(q,tmp3);
            multiply_bit(tmp2,tmp3,tmp);
            _base_n = BytesBase<T>(tmp);
            for (int i = 0;i<BYTESIZE; i++){
                if(tmp[i]){
                    _size = i;
                }
            }
            T bias[BYTESIZE]={0};
            bias[0]=1;
            minus_bit(tmp2,bias);
            minus_bit(tmp3,bias);
            int2bi(e, _value_e);
            T tmp1[BYTESIZE] = {0};
             multiply_bit(tmp2,tmp3,tmp1);
            _ok=inv(_value_e,tmp1,_key_d);
            if(!_ok){
                std::cout<<"init failed"<<std::endl;
            }
        }
        
        void encrp_data(T* text, T* result){
            _base_n.exp(text, _value_e, result);
        }
        void decrp_data(T* result, T* text){
            _base_n.exp(result, _key_d, text);
        }

        int encrp(T* a, int size, T* encode_text){
            T text[BYTESIZE]={0};
            T seg_result[BYTESIZE]={0};
            int encode_size = (sizeof(a)-1)/_size+1;
            int text_index = 0;
            int encode_index = 0;
            for(int i=0;i<size;i++){
                if(text_index==_size){
                    encrp_data(text, seg_result);
                    for (int j= 0; j <_size+1;j++){
                        encode_text[encode_index] = seg_result[j];
                        encode_index += 1;
                    }
                    text_index = 0;
                }
                text[text_index] = a[i];
                text_index += 1;
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

        void decrp(T *encode_text, int size, T* decode_text)
        {
            int decode_index = 0;
            int text_index = 0;
            T seg_encode[BYTESIZE]={0};
            T seg_result[BYTESIZE]={0};
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
        void decrp(char* o, int char_size, char* decode_text1 ){
            int base_size = char_size / sizeof(T) / 2;
            decode_text1[base_size / (_size + 1) * _size * sizeof(T)] = 0;
            T encode_test1[base_size];
            int size = str2char(o, char_size, encode_test1);
            decrp(encode_test1, size , (T*)decode_text1);
        }
        void encrp(const char *o, int char_size, char *result){
            int base_size = char_size / sizeof(T) + 1;
            T input[base_size];
            for (int i = 0; i < char_size; i++){
                *(((char *)input) + i) = o[i];
            }
            int decode_size = ((base_size-1)/_size+1)*(_size+1);
            T encode_text[decode_size];
            int size=encrp(input, base_size, encode_text);
            size = char2str(encode_text, size, result);
        }
        int get_size(){
            return _size;
        }

        int get_encode_size(int size){
            int base_size = size / sizeof(T) + 1;
            return  ((base_size-1)/_size+1)*(_size+1)*2*sizeof(T)+1;
        }

         int get_decode_size(int size){
            int base_size = size/sizeof(T)/2;
            return base_size / (_size + 1) * _size * sizeof(T) + 1;
         }
};




int main(){
    std::string a = "你好,你叫什么2121";
    Encrpyt<BASETYPE> enc(5,500887,500891);
    char result[enc.get_encode_size(a.size())];
    enc.encrp(a.c_str(), a.size(), result );
    std::cout<<result<<std::endl;
    
    char decode_text[enc.get_decode_size(sizeof(result))];
    enc.decrp(result,sizeof(result)-1, decode_text);
    std::cout << decode_text << std::endl;
    return 0;
}
