#!/bin/sh

sed 's/xa[khi]/material->Lookup[TILL_INDEX(irho,khi)].v/g' routine_old.txt > cmp
sed -i 's/xa[klo]/material->Lookup[TILL_INDEX(irho,klo)].v/g' cmp
sed -i 's/xa[k]/material->Lookup[TILL_INDEX(irho,k)].v/g' cmp

sed -i 's/ya[khi]/material->Lookup[TILL_INDEX(irho,khi)].u/g' cmp
sed -i 's/ya[klo]/material->Lookup[TILL_INDEX(irho,klo)].u/g' cmp

sed -i 's/y2a[khi]/material->Lookup[TILL_INDEX(irho,khi)].udv2/g' cmp
sed -i 's/y2a[klo]/material->Lookup[TILL_INDEX(irho,klo)].udv2/g' cmp

exit 1

sed 's/material->Lookup[TILL_INDEX(irho,khi)].v/xa[khi]/g' routine_new.txt > cmp
sed -i 's/material->Lookup[TILL_INDEX(irho,klo)].v/xa[klo]/g' cmp
sed -i 's/material->Lookup[TILL_INDEX(irho,k)].v/xa[k]/g' cmp

sed -i 's/material->Lookup[TILL_INDEX(irho,khi)].u/ya[khi]/g' cmp
sed -i 's/material->Lookup[TILL_INDEX(irho,klo)].u/ya[klo]/g' cmp

sed -i 's/material->Lookup[TILL_INDEX(irho,khi)].udv2/y2a[khi]/g' cmp
sed -i 's/material->Lookup[TILL_INDEX(irho,klo)].udv2/y2a[klo]/g' cmp

