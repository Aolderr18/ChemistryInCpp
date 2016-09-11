/*
 * The purpose of this program is to balance chemical equations.
 */

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
using namespace std;

/*
 * The following enumeration is used because the chemical elements are sequential
 * given that each next one has an atomic number exactly one higher than the 
 * previous one.
 */
enum element { None, Hydrogen, Helium, Lithium, Beryllium, Boron, Carbon, Nitrogen, Oxygen, Fluorine, Neon, Sodium, Magnesium, Aluminum, Silicon,
    Phosphorus, Sulfur, Chlorine, Argon, Potassium, Calcium, Scandium, Titanium, Vanadium, Chromium, Manganese, Iron, Cobalt, Nickel, Copper, Zinc, 
    Gallium, Germanium, Arsenic, Selenium, Bromine, Krypton, Rubidium, Strontium, Yttrium, Zirconium, Niobium, Molybdenum, Technetium, Ruthenium, Rhodium, 
    Palladium, Silver, Cadmium, Indium, Tin, Antimony, Tellurium, Iodine, Xenon, Cesium, Barium, Lanthanum, Cerium, Praseodymium, Neodymium, Promethium, 
    Samarium, Europium, Gadolinium, Terbium, Dysprosium, Holmium, Erbium, Thulium, Ytterbium, Lutetium, Hafnium, Tantalum, Tungsten, Rhenium, Osmium, 
    Iridium, Platinum, Gold, Mercury, Thallium, Lead, Bismuth, Polonium, Astatine, Radon, Francium, Radium, Actinium, Thorium, Protactinium, Uranium,
    Neptunium, Plutonium, Americium, Curium, Berkelium, Californium, Einsteinium, Fermium, Mendelevium, Nobelium, Lawrencium, Rutherfordium, Dubnium, 
    Seaborgium, Bohrium, Hassium, Meitnerium, Darmstadtium, Roentgenium, Ununbium, Ununtrium, Ununquadium, Ununpentium, Ununhexium, Ununseptium, Ununoctium };

map<element, string> elementNames(void)// This converts each element into a string.
{
    map<element, string> names;
    names[Hydrogen] = "Hydrogen";
    names[Helium] = "Helium";
    names[Lithium] = "Lithium";
    names[Beryllium] = "Beryllium";
    names[Boron] = "Boron";
    names[Carbon] = "Carbon";
    names[Nitrogen] = "Nitrogen";
    names[Oxygen] = "Oxygen";
    names[Fluorine] = "Fluorine";
    names[Neon] = "Neon";
    names[Sodium] = "Sodium";
    names[Magnesium] = "Magnesium";
    names[Aluminum] = "Aluminum";
    names[Silicon] = "Silicon";
    names[Phosphorus] = "Phosphorus";
    names[Sulfur] = "Sulfur";
    names[Chlorine] = "Chlorine";
    names[Argon] = "Argon";
    names[Potassium] = "Potassium";
    names[Calcium] = "Calcium";
    names[Scandium] = "Scandium";
    names[Titanium] = "Titanium";
    names[Vanadium] = "Vanadium";
    names[Chromium] = "Chromium";
    names[Manganese] = "Manganese";
    names[Iron] = "Iron";
    names[Cobalt] = "Cobalt";
    names[Nickel] = "Nickel";
    names[Copper] = "Copper";
    names[Zinc] = "Zinc";
    names[Gallium] = "Gallium";
    names[Germanium] = "Germanium";
    names[Arsenic] = "Arsenic";
    names[Selenium] = "Selenium";
    names[Bromine] = "Bromine";
    names[Krypton] = "Krypton";
    names[Rubidium] = "Rubidium";
    names[Strontium] = "Strontium";
    names[Yttrium] = "Yttrium";
    names[Zirconium] = "Zirconium";
    names[Niobium] = "Niobium";
    names[Molybdenum] = "Molybdenum";
    names[Technetium] = "Technetium";
    names[Ruthenium] = "Ruthenium";
    names[Rhodium] = "Rhodium";
    names[Palladium] = "Palladium";
    names[Silver] = "Silver";
    names[Cadmium] = "Cadmium";
    names[Indium] = "Indium";
    names[Tin] = "Tin";
    names[Antimony] = "Antimony";
    names[Tellurium] = "Tellurium";
    names[Iodine] = "Iodine";
    names[Xenon] = "Xenon";
    names[Cesium] = "Cesium";
    names[Barium] = "Barium";
    names[Lanthanum] = "Lanthanum";
    names[Cerium] = "Cerium";
    names[Praseodymium] = "Praseodymium";
    names[Neodymium] = "Neodymium";
    names[Promethium] = "Promethium";
    names[Samarium] = "Samarium";
    names[Europium] = "Europium";
    names[Gadolinium] = "Gadolinium";
    names[Terbium] = "Terbium";
    names[Dysprosium] = "Dysprosium";
    names[Holmium] = "Holmium";
    names[Erbium] = "Erbium";
    names[Thulium] = "Thulium";
    names[Ytterbium] = "Ytterbium";
    names[Lutetium] = "Lutetium";
    names[Hafnium] = "Hafnium";
    names[Tantalum] = "Tantalum";
    names[Tungsten] = "Tungsten";
    names[Rhenium] = "Rhenium";
    names[Osmium] = "Osmium";
    names[Iridium] = "Iridium";
    names[Platinum] = "Platinum";
    names[Gold] = "Gold";
    names[Mercury] = "Mercury";
    names[Thallium] = "Thallium";
    names[Lead] = "Lead";
    names[Bismuth] = "Bismuth";
    names[Polonium] = "Polonium";
    names[Astatine] = "Astatine";
    names[Radon] = "Radon";
    names[Francium] = "Francium";
    names[Radium] = "Radium";
    names[Actinium] = "Actinium";
    names[Thorium] = "Thorium";
    names[Protactinium] = "Protactinium";
    names[Uranium] = "Uranium";
    names[Neptunium] = "Neptunium";
    names[Plutonium] = "Plutonium";
    names[Americium] = "Americium";
    names[Curium] = "Curium";
    names[Berkelium] = "Berkelium";
    names[Californium] = "Californium";
    names[Einsteinium] = "Einsteinium";
    names[Fermium] = "Fermium";
    names[Mendelevium] = "Mendelevium";
    names[Nobelium] = "Nobelium";
    names[Lawrencium] = "Lawrencium";
    names[Rutherfordium] = "Rutherfordium";
    names[Dubnium] = "Dubnium";
    names[Seaborgium] = "Seaborgium";
    names[Bohrium] = "Bohrium";
    names[Hassium] = "Hassium";
    names[Meitnerium] = "Meitnerium";
    names[Darmstadtium] = "Darmstadtium";
    names[Roentgenium] = "Roentgenium";
    names[Ununbium] = "Uub";
    names[Ununtrium] = "Uut";
    names[Ununquadium] = "Uuq";
    names[Ununpentium] = "Uup";
    names[Ununhexium] = "Uuh";
    names[Ununseptium] = "Uus";
    names[Ununoctium] = "Uuo";
    return names;
}

/*
 * The following data structure links each element with its corresponding 
 * symbol.
 */
map<element, string> elementSymbols(void)
{
    map<element, string> symbols;
    symbols[Hydrogen] = "H";
    symbols[Helium] = "He";
    symbols[Lithium] = "Li";
    symbols[Beryllium] = "Be";
    symbols[Boron] = "B";
    symbols[Carbon] = "C";
    symbols[Nitrogen] = "N";
    symbols[Oxygen] = "O";
    symbols[Fluorine] = "F";
    symbols[Neon] = "Ne";
    symbols[Sodium] = "Na";
    symbols[Magnesium] = "Mg";
    symbols[Aluminum] = "Al";
    symbols[Silicon] = "Si";
    symbols[Sulfur] = "S";
    symbols[Chlorine] = "Cl";
    symbols[Argon] = "Ar";
    symbols[Potassium] = "K";
    symbols[Calcium] = "Ca";
    symbols[Scandium] = "Sc";
    symbols[Titanium] = "Ti";
    symbols[Vanadium] = "V";
    symbols[Chromium] = "Cr";
    symbols[Manganese] = "Mn";
    symbols[Iron] = "Fe";
    symbols[Cobalt] = "Co";
    symbols[Nickel] = "Ni";
    symbols[Copper] = "Cu";
    symbols[Zinc] = "Zn";
    symbols[Gallium] = "Ga";
    symbols[Germanium] = "Ge";
    symbols[Arsenic] = "As";
    symbols[Selenium] = "Se";
    symbols[Bromine] = "Br";
    symbols[Krypton] = "Kr";
    symbols[Rubidium] = "Rb";
    symbols[Strontium] = "Sr";
    symbols[Yttrium] = "Y";
    symbols[Zirconium] = "Zr";
    symbols[Niobium] = "Nb";
    symbols[Molybdenum] = "Mo";
    symbols[Technetium] = "Tc";
    symbols[Ruthenium] = "Ru";
    symbols[Rhodium] = "Rh";
    symbols[Palladium] = "Pd";
    symbols[Silver] = "Ag";
    symbols[Cadmium] = "Cd";
    symbols[Indium] = "In";
    symbols[Tin] = "Sn";
    symbols[Antimony] = "Sb";
    symbols[Tellurium] = "Te";
    symbols[Iodine] = "I";
    symbols[Xenon] = "Xe";
    symbols[Cesium] = "Cs";
    symbols[Barium] = "Ba";
    symbols[Lanthanum] = "La";
    symbols[Cerium] = "Ce";
    symbols[Praseodymium] = "Pr";
    symbols[Neodymium] = "Nd";
    symbols[Promethium] = "Pm";
    symbols[Samarium] = "Sm";
    symbols[Europium] = "Eu";
    symbols[Gadolinium] = "Gd";
    symbols[Terbium] = "Tb";
    symbols[Dysprosium] = "Dy";
    symbols[Holmium] = "Ho";
    symbols[Erbium] = "Er";
    symbols[Thulium] = "Tm";
    symbols[Ytterbium] = "Yb";
    symbols[Lutetium] = "Lu";
    symbols[Hafnium] = "Hf";
    symbols[Tantalum] = "Ta";
    symbols[Tungsten] = "W";
    symbols[Rhenium] = "Re";
    symbols[Osmium] = "Os";
    symbols[Iridium] = "Ir";
    symbols[Platinum] = "Pt";
    symbols[Gold] = "Au";
    symbols[Mercury] = "Hg";
    symbols[Thallium] = "Tl";
    symbols[Lead] = "Pb";
    symbols[Bismuth] = "Bi";
    symbols[Polonium] = "Po";
    symbols[Astatine] = "At";
    symbols[Radon] = "Rn";
    symbols[Francium] = "Fr";
    symbols[Radium] = "Ra";
    symbols[Actinium] = "Ac";
    symbols[Thorium] = "Th";
    symbols[Protactinium] = "Pa";
    symbols[Uranium] = "U";
    symbols[Neptunium] = "Np";
    symbols[Plutonium] = "Pu";
    symbols[Americium] = "Am";
    symbols[Curium] = "Cm";
    symbols[Berkelium] = "Bk";
    symbols[Californium] = "Cf";
    symbols[Einsteinium] = "Es";
    symbols[Fermium] = "Fm";
    symbols[Mendelevium] = "Md";
    symbols[Nobelium] = "No";
    symbols[Lawrencium] = "Lr";
    symbols[Rutherfordium] = "Rf";
    symbols[Dubnium] = "Db";
    symbols[Seaborgium] = "Sg";
    symbols[Bohrium] = "Bh";
    symbols[Hassium] = "Hs";
    symbols[Meitnerium] = "Mt";
    symbols[Darmstadtium] = "Ds";
    symbols[Roentgenium] = "Rg";
    symbols[Ununbium] = "Uub";
    symbols[Ununtrium] = "Uut";
    symbols[Ununquadium] = "Uuq";
    symbols[Ununpentium] = "Uup";
    symbols[Ununhexium] = "Uuh";
    symbols[Ununseptium] = "Uus";
    symbols[Ununoctium] = "Uuo";
    return symbols;
}

/*
 * The follow struct is to create instances of atoms, each with the element
 * symbol they represent and the amount of them. For instance O2 represents two
 * oxygen atoms.
 */
struct Atom {
    string elementSymbol;
    int count;
};

class Molecule {
private:
    vector<Atom> elements; /*
                            * These are all the atoms contained in a particular 
                            * instance of a molecule.
                            */
public:
    Molecule(void);
    Molecule(vector<Atom>);
    string userVisibleRepresentation(void); /*
                                             * This outputs the molecule's 
                                             * chemical formula.
                                             */
    vector<Atom> getElements(void);
    void setElements(vector<Atom>);
};

class ChemicalEquation {
private:
    vector<Molecule> leftSide;
    vector<Molecule> rightSide;
public:
    ChemicalEquation(void);
    ChemicalEquation(vector<Molecule> leftSide, vector<Molecule> rightSide);
    map<Molecule*, int> balance(void);
    vector<Molecule> getLeftSide(void);
    vector<Molecule> getRightSide(void);
    void setLeftSide(vector<Molecule>);
    void setRightSide(vector<Molecule>);
};

Molecule::Molecule(void)
{
    // Default Molecule constructor.
}

Molecule::Molecule(vector<Atom> elements)
{
    setElements(elements);
}

string Molecule::userVisibleRepresentation(void)
{
    stringstream chemicalFormula;
    for (int d = 0; d < elements.size(); ++d)
    {
        chemicalFormula << elements[d].elementSymbol;
        if (elements[d].count > 1)
            chemicalFormula << elements[d].count;
    }
    string userVisibleRepresentation = chemicalFormula.str();
    return userVisibleRepresentation;
}

vector<Atom> Molecule::getElements(void)
{
    return elements;
}

void Molecule::setElements(vector<Atom> elements)
{
    Molecule::elements = elements;
}

ChemicalEquation::ChemicalEquation(void)
{
    // Default chemical equation constructor
}

ChemicalEquation::ChemicalEquation(vector<Molecule> leftSide, vector<Molecule> rightSide)
{
    setLeftSide(leftSide);
    setRightSide(rightSide);
}

map<Molecule*, int> ChemicalEquation::balance(void)
{
    vector<string> allElementSymbols;
    bool contains;
    for (int i = 0; i < leftSide.size() + rightSide.size(); ++i)
        if (i < leftSide.size())
        {
            for (unsigned k = 0; k < leftSide[i].getElements().size(); ++k)
            {
                for (unsigned s = 0; s < allElementSymbols.size(); ++s)
                    if (allElementSymbols[s] == leftSide[i].getElements()[k].elementSymbol)
                        contains = true;
                if (!contains)
                    allElementSymbols.push_back(leftSide[i].getElements()[k].elementSymbol);
                contains = false;
            }
        }
        else
        {
            for (unsigned k = 0; k < rightSide[i - leftSide.size()].getElements().size(); ++k)
            {
                for (unsigned s = 0; s < allElementSymbols.size(); ++s)
                    if (allElementSymbols[s] == rightSide[i - leftSide.size()].getElements()[k].elementSymbol)
                        contains = true;
                if (!contains)
                    allElementSymbols.push_back(rightSide[i - leftSide.size()].getElements()[k].elementSymbol);
                contains = false;
            }
        }
    unsigned rows = 2 * allElementSymbols.size(), columns = leftSide.size() + rightSide.size() + 1;
    int** EQ_Matrix = new int * [rows];/*
                                        * This dynamic 2 dimensional array 
                                        * represents a linear algebra matrix
                                        * storing information about the molecules
                                        * in the equation.
                                        */
    unsigned n, h;
    for (n = 0; n < rows; ++n)
    {
        *(EQ_Matrix + n) = new int[columns];
        for (h = 0; h < columns; ++h)
            if (n % 2 == 0)
                *(*(EQ_Matrix + n) + h) = 0;
            else
                *(*(EQ_Matrix + n) + h) = 1;
    }
    for (h = 0; h < columns; ++h)
    {
        if (h < leftSide.size())
            for (int m = 0; m < leftSide[h].getElements().size(); ++m)
                for (int d = 0; d < allElementSymbols.size(); ++d)
                    if (leftSide[h].getElements()[m].elementSymbol == allElementSymbols[d])
                        *(*(EQ_Matrix + (d * 2)) + h) = leftSide[h].getElements()[m].count;
        if (h < columns - 1 and h >= leftSide.size()) 
        {
            for (int m = 0; m < rightSide[h - leftSide.size()].getElements().size(); ++m)
                for (int d = 0; d < allElementSymbols.size(); ++d)
                    if (rightSide[h - leftSide.size()].getElements()[m].elementSymbol == allElementSymbols[d])
                        *(*(EQ_Matrix + (d * 2)) + h) = rightSide[h - leftSide.size()].getElements()[m].count * -1;
        }
    }
    unsigned *pivotIndex = new unsigned[rows / 2];
    bool inEchelonForm = false;
    int x, y;
    while (!inEchelonForm) /*
                            * The row operations of the matrix should be performed
                            * until it is in echelon form.
                            */
    {
        for (x = 0; x < rows; x += 2)
        {
            for (y = 0; *(*(EQ_Matrix + x) + y) == 0 and y < columns; ++y);
            *(pivotIndex + (x / 2)) = y;
        }
        for (x = 2; x < rows; x += 2)
            if (pivotIndex[x / 2] < pivotIndex[(x - 2) / 2])
            {
                int* tempA = new int[columns];
                int* tempB = new int[columns];
                for (y = 0; y < columns; ++y)
                {
                    *(tempA + y) = *(*(EQ_Matrix + x) + y);
                    *(tempB + y) = *(*(EQ_Matrix + x + 1) + y);
                    *(*(EQ_Matrix + x) + y) = *(*(EQ_Matrix + x - 2) + y);
                    *(*(EQ_Matrix + x + 1) + y) = *(*(EQ_Matrix + x - 1) + y);
                    *(*(EQ_Matrix + x - 2) + y) = *(tempA + y);
                    *(*(EQ_Matrix + x - 1) + y) = *(tempB + y);
                }
                delete tempA;
                delete tempB;
            }
    inEchelonForm = true; /*
                           * The condition that determines if the loop should
                           * end should be set to 'true'. The following algorithm
                           * sets it back to 'false' if the conditions for being
                           * in echelon form are not met. If all of the conditions
                           * are met, the loop should end.
                           */
        for (x = 2; x < rows; x += 2)
            for (y = 0; y < x / 2; ++y)
                if (*(*(EQ_Matrix + x) + y) != 0)
                    inEchelonForm = false;
        if (inEchelonForm)
            continue;/*
                      * If the equation matrix is already in echelon form, there
                      * won't be any reason to perform any more row operations.
                      */
        int* rowOperationFraction = new int[2];
        for (x = 2; x < rows; ++x)
            for (y = 0; y < x / 2; ++y)
                if (*(*(EQ_Matrix + x) + y) != 0)
                {
                    *rowOperationFraction = *(*(EQ_Matrix + x) + y) * *(*(EQ_Matrix + x - 1) + y) * -1;
                    *(rowOperationFraction + 1) = *(*(EQ_Matrix + x + 1) + y) * *(*(EQ_Matrix + x - 2) + y);
                    int* tempA = new int[columns];
                    int* tempB = new int[columns];
                    for (int s = 0; s < columns; ++s)
                    {
                        *(tempA + s) = *(*(EQ_Matrix + x - 2) + s) * *rowOperationFraction;
                        *(tempB + s) = *(*(EQ_Matrix + x - 1) + s) * *(rowOperationFraction + 1);
                        int commonDenominator = *(*(EQ_Matrix + x + 1) + s) * *(tempB + s);
                        *(*(EQ_Matrix + x) + s) = *(*(EQ_Matrix + x) + s) * *(tempB + s) + *(*(EQ_Matrix + x + 1) + s) * *(tempA + s);
                        *(*(EQ_Matrix + x + 1) + s) = commonDenominator;
                    }
                    delete tempA;
                    delete tempB;
                    x = rows;
                    y = x;
                }
        delete rowOperationFraction;
    }
    int nonFreeVariables = rows / 2;
    int freeVariables = columns - (rows / 2) - 1;
    int** variableSolutionsNumer = new int * [rows / 2];
    int** variableSolutionsDenom = new int * [rows / 2];
    for (x = 0; x < rows; x += 2)
    {
        *(variableSolutionsNumer + (x / 2)) = new int[nonFreeVariables - (x / 2)];
        *(variableSolutionsDenom + (x / 2)) = new int[nonFreeVariables - (x / 2)];
        int* ratio = new int[2];
        *ratio = *(*(EQ_Matrix + x + 1) + (x / 2));
        *(ratio + 1) = *(*(EQ_Matrix + x) + (x / 2));
        for (int s = 0; s < (nonFreeVariables - (x / 2)); ++s)
        {
            *(*(variableSolutionsNumer + (x / 2)) + s) = *ratio * *(*(EQ_Matrix + x) + ((x / 2) + s + 1)) * -1;
            *(*(variableSolutionsDenom + (x / 2)) + s) = *(ratio + 1) * *(*(EQ_Matrix + x + 1) + ((x / 2) + s + 1));
        }
        delete ratio;
    }
    int** solutionNumeratorSet = new int * [rows / 2]; 
    int** solutionDenominatorSet = new int * [rows / 2];/*
                                                         * The algorithm used to compute the variable solutions
                                                         * converts the integer matrix elements into fractions.
                                                         */
    for (x = (rows / 2) - 1; x + 1; --x)
    {
        int o = (rows / 2) - 1 - x;
        *(solutionNumeratorSet + x) = new int[(rows / 2) - x];
        for (y = 0; y < (rows / 2) - x - freeVariables; ++y)
        {
            --o;
            for (int e = 0; e < freeVariables; ++e)
            {
                int tempNum = *(*(variableSolutionsNumer + x) + y) * *(*(variableSolutionsNumer + x + y + 1) + o + e);
                int tempDen = *(*(variableSolutionsDenom + x) + y) * *(*(variableSolutionsDenom + x + y + 1) + o + e);
                int commonDenominator = *(*(variableSolutionsDenom + x) + (rows / 2) - x - freeVariables) * tempDen;
                *(*(variableSolutionsNumer + x) + (rows / 2) - x - freeVariables) *= tempDen;
                *(*(variableSolutionsNumer + x) + (rows / 2) - x - freeVariables) += *(*(variableSolutionsDenom + x) + (rows / 2) - x - freeVariables) * tempNum;
                *(*(variableSolutionsDenom + x) + (rows / 2) - x - freeVariables) = commonDenominator;
            }
        }
    }
    for (x = 0; x < rows / 2; ++x)
        for (y = 0; y < freeVariables; ++y)
        {
            int numer = *(*(variableSolutionsNumer + x) + y + ((rows / 2) - 1 - x));
            int denom = *(*(variableSolutionsDenom + x) + y + ((rows / 2) - 1 - x));
            int divisor = 2;
            int destination;
            if (numer > denom)
                destination = denom;
            else
                destination = numer;
            while (divisor < destination)
            {
                if (numer % divisor == 0 and denom % divisor == 0)
                {
                    numer /= divisor;
                    denom /= divisor;
                    --divisor;
                }
                ++divisor;
            }
            *(*(variableSolutionsNumer + x) + y + ((rows / 2) - 1) - x) = numer;
            *(*(variableSolutionsDenom + x) + y + ((rows / 2) - 1) - x) = denom;
        }
    int max = 0;
    int testAmount = 0;
    bool lcmFound = false;
    for (int u = 0; !lcmFound; ++u)
        for (int w = 0; w < freeVariables; ++w)
        {
            for (x = 0; x < rows / 2; ++x)
            {
                if (u == 0)
                {
                    if (*(*(variableSolutionsDenom + x) + (rows / 2) - 1 - x) > max)
                    {
                        max = *(*(variableSolutionsDenom + x) + nonFreeVariables - 1 + w);
                        testAmount = max;
                    }
                    continue;
                }
                if (testAmount % *(*(variableSolutionsDenom + x) + (rows / 2) - 1 - x) != 0)
                {
                    testAmount += max;
                    lcmFound = false;
                    continue;
                }
                lcmFound = true;
            }
        }
    vector<int> amounts;/*
                         * This integer vector represents the amount of each 
                         * molecule in the chemical equation should be used to
                         * properly balance the equation.
                         */
    for (x = 0; x < rows / 2; ++x)
        amounts.push_back((*(*(variableSolutionsNumer + x) + nonFreeVariables - x - 1) * testAmount) / *(*(variableSolutionsDenom + x) + nonFreeVariables - x - 1));
    for (int d = 0; d < freeVariables; ++d)
        amounts.push_back(testAmount);
    delete variableSolutionsNumer;
    delete variableSolutionsDenom;
    map<Molecule*, int> bal;
    for (int v = 0; v < amounts.size(); ++v)
        if (v < leftSide.size())
            bal.insert(std::make_pair(new Molecule(leftSide[v]), amounts[v]));
        else
            bal.insert(std::make_pair(new Molecule(rightSide[v - leftSide.size()]), amounts[v]));
    return bal;
}

vector<Molecule> ChemicalEquation::getLeftSide(void)
{
    return leftSide;
}

vector<Molecule> ChemicalEquation::getRightSide(void)
{
    return rightSide;
}

void ChemicalEquation::setLeftSide(vector<Molecule> leftSide)
{
    ChemicalEquation::leftSide = leftSide;
}

void ChemicalEquation::setRightSide(vector<Molecule> rightSide)
{
    ChemicalEquation::rightSide = rightSide;
}

int main(void)
{
    ofstream outFile;
    outFile.open("ffile.dat", ios::out);/*
                                           * This is the file that the information about
                                           * the balanced chemical equation will be 
                                           * written to.
                                           */
    Atom firstCarbon, firstHydrogen, secondOxygen, thirdCarbon, thirdOxygen, fourthHydrogen, fourthOxygen;
    firstCarbon.elementSymbol = "C";
    firstCarbon.count = 3;
    firstHydrogen.elementSymbol = "H";
    firstHydrogen.count = 8;
    secondOxygen.elementSymbol = "O";
    secondOxygen.count = 2;
    thirdCarbon.elementSymbol = "C";
    thirdCarbon.count = 1;
    thirdOxygen.elementSymbol = "O";
    thirdOxygen.count = 2;
    fourthHydrogen.elementSymbol = "H";
    fourthHydrogen.count = 2;
    fourthOxygen.elementSymbol = "O";
    fourthOxygen.count = 1;
    vector<Atom> farLeftAtoms, nearLeftAtoms, nearRightAtoms, farRightAtoms;
    farLeftAtoms.push_back(firstCarbon);
    farLeftAtoms.push_back(firstHydrogen);
    nearLeftAtoms.push_back(secondOxygen);
    nearRightAtoms.push_back(thirdCarbon);
    nearRightAtoms.push_back(thirdOxygen);
    farRightAtoms.push_back(fourthHydrogen);
    farRightAtoms.push_back(fourthOxygen);
    Molecule farLeftMolecule, nearLeftMolecule, nearRightMolecule, farRightMolecule;
    farLeftMolecule.setElements(farLeftAtoms);
    nearLeftMolecule.setElements(nearLeftAtoms);
    nearRightMolecule.setElements(nearRightAtoms);
    farRightMolecule.setElements(farRightAtoms);
    vector<Molecule> leftSideOfEquation, rightSideOfEquation;
    leftSideOfEquation.push_back(farLeftMolecule);
    leftSideOfEquation.push_back(nearLeftMolecule);
    rightSideOfEquation.push_back(nearRightMolecule);
    rightSideOfEquation.push_back(farRightMolecule);
    ChemicalEquation chemEQ;
    chemEQ.setLeftSide(leftSideOfEquation);
    chemEQ.setRightSide(rightSideOfEquation);
    /*
     * In this example, the chemical equation that will be balanced will use 
     * four molecules, namely propane(C3H8), oxygen(O2), carbon dioxide(CO2),
     * and water(H2O).
     */
    map<Molecule*, int> howManyOfEach = chemEQ.balance();
                                      /*
                                       * This mapping data structure links each
                                       * molecule in the chemical equation with
                                       * the lowest integer count they should be
                                       * assigned to properly balance the equation.
                                       */
    map<element, string> names = elementNames();
    outFile << "To properly balance the chemical equation, the coefficients ";
    for (map<Molecule*, int>::iterator it = howManyOfEach.begin(); it != howManyOfEach.end(); ++it)
    {
        outFile << it->second << " for " << it->first->userVisibleRepresentation() << ",\n";
    }
    outFile << " must be used.";
    outFile.close();
    ifstream inFile; /*
                      * The cout stream will read data from the previously 
                      * file with data output to it.
                      */
    inFile.open("ffile.dat");
    char data[100];
    bool lastWord = false;/*
                           * The file will continue to be read from until all the
                           * data has been read. The final word must be "used".
                           */
    do
    {
        inFile >> data;
        cout << data << " ";
        stringstream word;
        for (int p = 0; p < 5; ++p)
            word << data[p];
        if (word.str() == "used.")
            lastWord = true;
    }
    while (!lastWord);
    inFile.close();
    return 0;
}
