#include "Amino.h"


Amino::Amino(Direction dir, bool isHydrophob)
{
	this->Dir = dir;
	this->Hydrophob = isHydrophob;
}

Amino::Amino(void)
{
	this->Dir = FORWARD;
}

Amino::~Amino(void)
{
}

// GETTER
Direction Amino::GetDir()
{
	return this->Dir;
}

bool Amino::IsHydrophob()
{
	return this->Hydrophob;
}

// SETTER
void Amino::SetDir(Direction dir)
{
	this->Dir = dir;
}

void Amino::SetHydrophob(bool hydrophob)
{
	this->Hydrophob = hydrophob;
}