#pragma once
#define BASETYPE unsigned short
#define BYTESIZE 128/sizeof(BASETYPE)/8
template<typename T>
void print(T* a){
    for (int i = 0; i < BYTESIZE;i++){
        std::cout << a[i]-0 << " ";
    }
    std::cout << std::endl;
}

template<typename T>
bool compare(const T* data, const T* bytes){
    for (int i = BYTESIZE-1; i >= 0; i--){
        if(data[i]>bytes[i]){
            return true;
        }else if(data[i]<bytes[i]){
            return false;
        }
    }
    return false;
}

template<typename T>
bool equal(const T* data, const T* bytes){
    for (int i = 0; i < BYTESIZE; i++){
        if(data[i]!=bytes[i]){
            return false;
        }
    }
    return true;
}

template<typename T>
bool equal_zero(const T* data){
    for (int i = 0; i < BYTESIZE; i++){
        if(data[i]!=0){
            return false;
        }
    }
    return true;
}

template<typename T>
bool char_self_add(T& a, T b, bool add){
  a+=b;
  if(add){
     a+=1;
    return a<=b;
   }
    else{
	return a<b; 
   }  
}

template<typename T>
bool har_self_add(T& a, T b, bool add){
    bool overflow = false;
   T tmp = 0;
    while(b){
        tmp = b;
        b &= a;
        if(b >>(sizeof(T)*8-1)){
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
            if(b >>(sizeof(T)*8-1)){
                overflow = true;
            }
            b <<=  1;
            a ^= tmp;
        }
    }
    return overflow;
}

template<typename T>
bool add_bit(T* a,T *b){
    bool add = false;
    for (int i = 0; i < BYTESIZE; i++){
        add = char_self_add(a[i], b[i], add);
    }
    return add;
}
template<typename T>
bool minus_bit(T* a,T *b){
    bool add = true;
    for (int i = 0; i < BYTESIZE; i++){
        add = char_self_add(a[i], (T)~b[i], add);
    }
    return add;
}

template<typename T>
int log2(T* a){
     int size = 0;
    int start = BYTESIZE*8*sizeof(T)-1;
    int first =0;
    int second =0;
    while (size==0&&start>=0){
        first=start&(8*sizeof(T)-1);
        second = start/(sizeof(T)*8);
        if(a[second]&(((T)1)<<first)){
            return start;
        }
        start-=1;
    }
    return -1;
}

template<typename T>
void left_bit(T* a,int size){
    int pre_fist=(BYTESIZE-1-(size/(sizeof(T)*8)));
    int pre_bias =(size&(8*sizeof(T)-1));
    T pre=0;
    for(int i=BYTESIZE-1;i>=0;i--){
        if(pre_fist<0){
            a[i]=0;
        }else if(pre_fist==0){
            a[i]=(a[pre_fist]<<pre_bias);
        }else{
            a[i]=(a[pre_fist]<<pre_bias)+(a[pre_fist-1]>>(8*sizeof(T)-pre_bias));
        }
        pre_fist-=1;
    }
}

template<typename T>
void right_bit(T* a,int size){
    int pre_fist=(size/(sizeof(T)*8));
    int pre_bias =(size&(8*sizeof(T)-1));
    T pre=0;
    for(int i=0;i<BYTESIZE;i++){
        if(pre_fist>BYTESIZE){
            a[i]=0;
        }else if(pre_fist==BYTESIZE-1){
            a[i]=(a[pre_fist]>>(pre_bias));
        }else{
            a[i]=(a[pre_fist]>>(pre_bias))+(a[pre_fist+1]<<(8*sizeof(T)-pre_bias));
        }
        pre_fist+=1;
    }
}

template<typename T>
bool divide_bit(T* a,T *b,T *redis){
    T ans[BYTESIZE] ={0};
    T bias[BYTESIZE];
    for(int i=0;i<BYTESIZE;i++){
        bias[i]=b[i];
    }
    int size = log2(a)-log2(b);
    if(size<0){
        for(int i=0;i<BYTESIZE;i++){
            redis[i]=a[i];
            a[i]=0;
        }
        return true;
    }
    left_bit(bias,size);

    while(!compare(b,a)){
        
        if(!compare(bias,a)){
            ans[size/(sizeof(T)*8)]+=(((T)1)<<(size&(8*sizeof(T)-1)));
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

template<typename T>
bool 
multiply_bit(T* a,T *b, T*ans){
    bool a_flat[BYTESIZE*8*sizeof(T)];
    bool b_flat[BYTESIZE*8*sizeof(T)];
    for(T i=0;i<BYTESIZE*8*sizeof(T);i++){
        T first =i/(sizeof(T)*8);
        T second =i&(8*sizeof(T)-1);
        if((a[first]>>second)&1){
            a_flat[i]=true;
        }else{
            a_flat[i]=false;
        }
        if((b[first]>>second)&1){
            b_flat[i]=true;
        }else{
            b_flat[i]=false;
        }

    }
    T  first=0;
    T second=0;
    T number = 1;
    T bias[BYTESIZE]={0};
    T tmpve[BYTESIZE]={0};
    tmpve[BYTESIZE-1]=12;
    while(first<BYTESIZE){
        bias[first]=(((T)1)<<second);
        for(T i=0;i<number;i++){
            if(a_flat[i]&&b_flat[number-1-i]){
                if(add_bit(ans, bias)){
                    return true;
                }
            }
        }
        second+=1;
        number+=1;
        if(second==8*sizeof(T)){
            bias[first]=0;
            second = 0;
            first+=1;
            if(first<BYTESIZE){
                bias[first]=1;
            }
        }
    }
    T start =1;
    while (number<2*BYTESIZE*8*sizeof(T)) {
        for(T i=start;i<number;i++){
            if(a_flat[i]&&b_flat[number-1-i]){
                return true;
            }
        }
        number+=1;
    }
    
    return false;
}

template<typename T>
bool ex_gcd(T* a, T* b, T* x, T* y, bool& x_flag, bool& y_flag){
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
        T next_a[BYTESIZE];
        T next_b[BYTESIZE];
        T resi[BYTESIZE];
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

template<typename T>
bool inv(T* t, T* p, T* ans){
    T x[BYTESIZE];
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
