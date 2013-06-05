#pragma once
#include <iostream>

class Amino;

class Protein
{
	private:
		Amino* Aminos;
		int Length;

		// Caching Attribute
		int Fitness;
		std::string Id;

	public:
		Protein(int length);
		Protein(std::string);
		Protein();
		~Protein();

		Amino* GetAminos();

		int CalcFitness();
		Protein* Copy();
		std::string Identifier();
		void Mutate(double mutationRate);
		void Show();
};