//
// Created by juzix on 2021/6/1.
//

#ifndef HYRAX_P224_TYPEDEF_HPP
#define HYRAX_P224_TYPEDEF_HPP

typedef unsigned __int128 u128;
typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned char u8;

typedef __int128 i128;
typedef long long i64;
typedef int i32;
typedef char i8;

typedef u64 p224_limb;
typedef u128 p224_widelimb;

typedef p224_limb p224_felem[4];
typedef p224_widelimb p224_widefelem[7];

#endif //HYRAX_P224_TYPEDEF_HPP
