Sets
    year_all bla bla /2010, 2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050,2008/
;

Alias(year_all,vintage);
Alias(year_all,year_all2);
Alias(year_all,year_all3);

Sets
         seq_period(year_all,year_all2)    mapping of one period ('year_all') to the next ('year_all2')

Parameter

    year_order(year_all)       order for members of set 'year_all'
;

seq_period(year_all,year_all2)$( ORD(year_all) + 1 = ORD(year_all2) ) = yes ;
year_order(year_all) = ORD(year_all) ;

display year_order
display seq_period