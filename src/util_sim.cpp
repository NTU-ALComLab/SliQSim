#include "util_sim.h"


/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1(int length, int *reg)
{
    int one = 1, carry = 0, i;
    for (i = 0; i < length; i++)
    {
        carry = reg[i] & one;
        reg[i] = reg[i] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1_start(int length, int *reg, int start)
{
    int one = 1, carry = 0, i;
    for (i = start; i < length; i++)
    {
        carry = reg[i] & one;
        reg[i] = reg[i] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void full_adder_plus_1_measure(int length, int *reg, int *order)
{
    int one = 1, carry = 0, i;
    for (i = length - 1; i >= 0; i--)
    {
        carry = reg[order[i]] & one;
        reg[order[i]] = reg[order[i]] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int int_array_full_check(int length, int *reg)
{
    int check = 1, i;
    for (i = 0; i < length; i++)
        check *= reg[i];

    return check;
}