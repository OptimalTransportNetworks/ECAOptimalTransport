
use "/Users/sebastiankrantz/Documents/ECAOptimalTransport/data/ECA_database_fullsample.dta", clear
* su dsalesusd dvausd vaprodusd dwagebillusd dcapusd denergyusd dmatusd emp
keep if !missing(dsalesusd, vaprodusd, emp) & (!missing(dwagebillusd) | !missing(dcapusd) | !missing(denergyusd) | !missing(dmatusd))

gen n_firms = 1
label variable n_firms  "Number of Firms"

foreach var in dsalesusd vaprodusd emp dwagebillusd dcapusd denergyusd dmatusd {
	gen n_zero_`var' = missing(`var')
}

* Save labels before collapse
foreach v of var * {
  local l`v' : variable label `v'
}

collapse (sum) dsalesusd dvausd vaprodusd dwagebillusd dcapusd denergyusd dmatusd emp exporter exports imports imports_value n_firms-n_zero_dmatusd, by(countrycode countryname country year shapeISO nace1d)

* Restore labels
foreach v of var * {
  label variable `v' `"`l`v''"'
}

gen single_firm = n_firms == 1
table single_firm

* Rsh=rev.share; N=nuts; T=total;   
* Revenue shares (for each input)=expenditures/sales 
 
* Total sales, value-added, wagebill, capital, energy expenses, material expenses
bysort country shapeISO year: egen NTSales=sum(dsalesusd) 
bysort country shapeISO year: egen NTVa=sum(vaprodusd)
bysort country shapeISO year: egen NTTCost=sum(dwagebillusd + dcapusd + denergyusd + dmatusd)
bysort country shapeISO year: egen NTW=sum(dwagebillusd) 
bysort country shapeISO year: egen NTFA=sum(dcapusd)
bysort country shapeISO year: egen NTEn=sum(denergyusd)
bysort country shapeISO year: egen NTRMat=sum(dmatusd)
bysort country shapeISO year: egen NTL=sum(emp)

label variable NTSales      "Total Sales (US$, PPP-adjusted) at NUTS level"
label variable NTVa         "Total Value-Added (US$, PPP-adjusted) at NUTS level"
label variable NTTCost      "Total Cost of Production (Wages + Capital + Energy + Materials; US$, PPP-adjusted) at NUTS level"
label variable NTW          "Total Wage Bill (US$, PPP-adjusted) at NUTS level"
label variable NTFA         "Total Fixed Assets / Capital Costs (US$, PPP-adjusted) at NUTS level"
label variable NTEn         "Total Energy Expenses (US$, PPP-adjusted) at NUTS level"
label variable NTRMat       "Total Raw Material Expenses (US$, PPP-adjusted) at NUTS level"
label variable NTL          "Total Employment (Number of Employees) at NUTS level"

* Total Shares
gen RshW=NTW/NTSales /* Wagebill (US$; PPP-adj),  Sales (US$; PPP-adj) */
gen RshFA=NTFA/NTSales /* Fixed assets (US$; PPP-adj) */
gen RshEn=NTEn/NTSales /* Energy (US$; PPP-adj) */
gen RshMat=NTRMat/NTSales  /* Raw input cost (US$, PPP adjusted) */
gen RshTCost=NTTCost/NTSales /* Total cost of production (US$, PPP adjusted) */

label variable RshW         "Share of Wage Bill in Total Sales (Wages/Sales, US$, PPP-adjusted) at NUTS level"
label variable RshFA        "Share of Fixed Assets in Total Sales (Capital/Sales, US$, PPP-adjusted) at NUTS level"
label variable RshEn        "Share of Energy Expenses in Total Sales (Energy/Sales, US$, PPP-adjusted) at NUTS level"
label variable RshMat       "Share of Material Expenses in Total Sales (Materials/Sales, US$, PPP-adjusted) at NUTS level"
label variable RshTCost     "Share of Total Production Costs in Total Sales (Total Costs/Sales, US$, PPP-adjusted) at NUTS level"

* Sectorial variables (cost shares = production function elasticities) at NUTS level => Likely not needed
gen NTTCost_nace1d=dwagebillusd + dcapusd + denergyusd + dmatusd
gen NshTCost_nace1d_C=NTTCost_nace1d/NTTCost
gen NshTCost_nace1d_S=NTTCost_nace1d/NTSales

label variable NTTCost_nace1d      "Total Sectoral Cost of Production (Wages + Capital + Energy + Materials; US$, PPP-adjusted) for NACE 1-digit sectors at NUTS level"
label variable NshTCost_nace1d_C   "Share of Sectoral Cost in Total Regional Production Cost (Sector Cost/Total Cost, US$, PPP-adjusted) at NUTS level"
label variable NshTCost_nace1d_S   "Share of Sectoral Cost in Total Regional Sales (Sector Cost/Sales, US$, PPP-adjusted) at NUTS level"

foreach var in dwagebillusd dcapusd denergyusd dmatusd {
	gen N_`var'_TCsh_nace1d = `var' / NTTCost_nace1d
}

* Sectorial variables (cost shares = production function elasticities) at Country level
bysort country nace1d year: egen CTTCost_nace1d=sum(dwagebillusd + dcapusd + denergyusd + dmatusd)
bysort country year: egen CTTCost=sum(dwagebillusd + dcapusd + denergyusd + dmatusd)
bysort country year: egen CTSales=sum(dsalesusd) 
gen CshTCost_nace1d_C=CTTCost_nace1d/CTTCost
gen CshTCost_nace1d_S=CTTCost_nace1d/CTSales

label variable CTTCost      "Total Cost of Production (Wages + Capital + Energy + Materials; US$, PPP-adjusted) at Country level"
label variable CTSales      "Total Sales (US$, PPP-adjusted) at Country level"
label variable CTTCost_nace1d      "Total Sectoral Cost of Production (Wages + Capital + Energy + Materials; US$, PPP-adjusted) for NACE 1-digit sectors at Country level"
label variable CshTCost_nace1d_C   "Share of Sectoral Cost in Total Regional Production Cost (Sector Cost/Total Cost, US$, PPP-adjusted) at Country level"
label variable CshTCost_nace1d_S   "Share of Sectoral Cost in Total Regional Sales (Sector Cost/Sales, US$, PPP-adjusted) at Country level"

foreach var in dwagebillusd dcapusd denergyusd dmatusd {
	bysort country nace1d year: egen C_`var'_nace1d=sum(`var')
	gen C_`var'_TCsh_nace1d = C_`var'_nace1d / CTTCost_nace1d
	drop C_`var'_nace1d
}

* Sectorial variables (cost shares) at Aggregate (Overall) Level
gen dtotalcost = dwagebillusd + dcapusd + denergyusd + dmatusd
gen dtotalcost_nzero = dtotalcost if n_firms != n_zero_dwagebillusd
gen dsalesusd_nzero = dsalesusd if n_firms != n_zero_dwagebillusd
bysort nace1d year: egen OTTCost_nace1d = sum(dtotalcost_nzero)
bysort year: egen OTTCost = sum(dtotalcost_nzero)
bysort year: egen OTSales = sum(dsalesusd_nzero) 
gen OshTCost_nace1d_C=OTTCost_nace1d/OTTCost
gen OshTCost_nace1d_S=OTTCost_nace1d/OTSales
drop dtotalcost_nzero dsalesusd_nzero

label variable OTTCost      "Total Cost of Production (Wages + Capital + Energy + Materials; US$, PPP-adjusted)"
label variable OTSales      "Total Sales (US$, PPP-adjusted) at Country level"
label variable OTTCost_nace1d      "Total Sectoral Cost of Production (Wages + Capital + Energy + Materials; US$, PPP-adjusted) for NACE 1-digit sectors"
label variable OshTCost_nace1d_C   "Share of Sectoral Cost in Total Regional Production Cost (Sector Cost/Total Cost, US$, PPP-adjusted)"
label variable OshTCost_nace1d_S   "Share of Sectoral Cost in Total Regional Sales (Sector Cost/Sales, US$, PPP-adjusted)"

foreach var in dwagebillusd dcapusd denergyusd dmatusd {
	gen `var'_nzero = `var' if n_firms != n_zero_dwagebillusd
	bysort nace1d year: egen O_`var'_nace1d=sum(`var'_nzero)
	gen O_`var'_TCsh_nace1d = O_`var'_nace1d / OTTCost_nace1d
	drop O_`var'_nace1d `var'_nzero
}

/* 
We then take these industry specific cost shares and aggregate them to a NUTS level by taking the sum of 
each industry cost share times that industry's share of revenue in the NUTS. This final calculation should 
deliver NUTS level cost shares that reflect the industry variation within that location.
*/

* Cost shares computed at country level
bysort country shapeISO year: egen CshTCost_nace1d_C_N = sum(CshTCost_nace1d_C * (dsalesusd / NTSales))
bysort country shapeISO year: egen C_dwagebillusd_TCsh_nace1d_N = sum(C_dwagebillusd_TCsh_nace1d * (dsalesusd / NTSales))
bysort country shapeISO year: egen C_dcapusd_TCsh_nace1d_N = sum(C_dcapusd_TCsh_nace1d * (dsalesusd / NTSales))
bysort country shapeISO year: egen C_denergyusd_TCsh_nace1d_N = sum(C_denergyusd_TCsh_nace1d * (dsalesusd / NTSales))
bysort country shapeISO year: egen C_dmatusd_TCsh_nace1d_N = sum(C_dmatusd_TCsh_nace1d * (dsalesusd / NTSales))

* Cost shares computed overall
bysort country shapeISO year: egen OshTCost_nace1d_C_N = sum(OshTCost_nace1d_C * (dsalesusd / NTSales))
bysort country shapeISO year: egen O_dwagebillusd_TCsh_nace1d_N = sum(O_dwagebillusd_TCsh_nace1d * (dsalesusd / NTSales))
bysort country shapeISO year: egen O_dcapusd_TCsh_nace1d_N = sum(O_dcapusd_TCsh_nace1d * (dsalesusd / NTSales))
bysort country shapeISO year: egen O_denergyusd_TCsh_nace1d_N = sum(O_denergyusd_TCsh_nace1d * (dsalesusd / NTSales))
bysort country shapeISO year: egen O_dmatusd_TCsh_nace1d_N = sum(O_dmatusd_TCsh_nace1d * (dsalesusd / NTSales))

* Optionally: collapse years as in paper (set to 1)
if 0 {
	foreach v of var * {
	  local l`v' : variable label `v'
	}
	collapse (mean) dsalesusd-O_dmatusd_TCsh_nace1d_N, by(countrycode countryname country shapeISO nace1d)
	foreach v of var * {
	  label variable `v' `"`l`v''"'
	}
}

* Wedges and gaps should be
foreach v in C O {
	gen gap_`v'_TCost = `v'shTCost_nace1d_C_N - RshTCost
	gen wedge_`v'_TCost = gap_`v'_TCost / RshTCost
	
	gen gap_`v'_W = `v'_dwagebillusd_TCsh_nace1d_N - RshW
	gen wedge_`v'_W = gap_`v'_W / RshW

	gen gap_`v'_FA = `v'_dcapusd_TCsh_nace1d_N - RshFA
	gen wedge_`v'_FA = gap_`v'_FA / RshFA

	gen gap_`v'_En = `v'_denergyusd_TCsh_nace1d_N - RshEn
	gen wedge_`v'_En = gap_`v'_En / RshEn

	gen gap_`v'_Mat = `v'_dmatusd_TCsh_nace1d_N - RshMat
	gen wedge_`v'_Mat = gap_`v'_Mat / RshMat
}

save "/Users/sebastiankrantz/Documents/ECAOptimalTransport/data/ECA_database_shapeISO.dta", replace
