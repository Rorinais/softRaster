#pragma once
#include<Windows.h>
#include<iostream>
#include<math.h>
#include<vector>
#include<array>
#include <cstdint>
#pragma comment(linker,"/subsystem:console /entry:wWinMainCRTStartup")
#define PI 3.1415

struct RGBA {
	byte mB;
	byte mG;
	byte mR;
	byte mA;
	RGBA(byte r = 255, byte g = 255, byte b = 255, byte a = 255){
		mR = r;
		mG = g;
		mB = b;
		mA = a;
	}
};