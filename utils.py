import pandas as pd 
import numpy as np 
import os
from scipy import optimize
import matplotlib.pyplot as plt

def kelvin(Tc): 
    """Convert degrees Celcius to degrees Kelvin
    Args: 
        - Tc (float): temperature in Celcius
    Returns: 
        - Tk (float): temperature in Kelvin
    """
    Tk = Tc + 273.15
    return Tk


def Lambda(Tc): 
    """Calculate latent heat of vaporization
    Args: 
        - Tc (float): Temperature in Celcius
    Returns: 
        - Lv (float): Latent heat of vaporization at input temperature in units of J/kg 
    """
    Lv = 3149000 - 2370*kelvin(Tc) #J/kg   
    return Lv


def satVapor(Tc):
    """Calculate saturation vapor pressure at a given temperature
    Args: 
        - Tc (float): Temperature in Celcius
    Returns: 
        - Vp (float): Saturation vapor pressure in units of Pa
    """
    A = 611 #Pa
    B = -5310 #K
    T0 = 273 #K
    Tk = kelvin(Tc) #K
    Vp = A*np.exp(B*((1/Tk)-(1/T0))) #Pa
    return Vp


def slope(Tc): 
    """Calculate the slope of the sat vapor curve 
    Args:
        - Tc (float): Temperature in Celcius
    Returns: 
        - s (float): Slope at Tc in units of Pa/K
    """
    C = 5310
    T = kelvin(Tc) #K
    s = (C*satVapor(Tc))/T**2 #Pa/K
    return s


def gamma(P,Tc):
    """Calculate the psychronmetric constant
    Args: 
        - Tc (float): Temperature in Celcius
        - P (float): pressure in kPa
    Returns: 
        - gamma (float): psychronmetric constant in units of Pa/K
    """
    Cp = 1005 #molar specific heat of dry air, J/(kg*K)
    Pa = 1000*P #convert kPa to Pa
    gamma = (Cp*Pa)/(Lambda(Tc)*0.622) #Pa/K
    return gamma


def E_pot(Tc, P, NETRAD, G): 
    """
    Args: 
        - Tc (float): Temperature in Celcius
        - P (float): Pressure in kPa
        - NETRAD (float): Net radiation in W/m^2
        - G (float): Soil heat flux in W/m^2
    Returns: 
        - E_pot (float): potential evapotranspiration in units of W/m^2
    """
    s = slope(Tc)
    p = gamma(P,Tc)
    E_pot = 1.26*(s/(s+p))*(NETRAD - G) #W/m^2
    return E_pot


def addESI_PT(fluxMonthly, upperLimit = 1.1): 
    """Calculate ESI PT and add as a column to the table 
    Args: 
        - fluxMonthly (pandas DataFrame): dataframe containing fluxnet monthly data
        - upperLimit (int, optional): upper limit to restrict ESI ratio (default to 1.1)
    Returns: 
        - fluxDF (pandas DataFrame): dataframe with new 'ESI PT (W/m^2)' column
    """
    fluxDF = fluxMonthly.copy()

    #pull data from columns in table
    Tc = fluxDF['TA (deg C)'].values
    P = fluxDF['PA (kPa)'].values
    NETRAD = fluxDF['NETRAD (W/m^2)'].values
    G = fluxDF['G (W/m^2)'].values
    LE = fluxDF['LE (W/m^2)'].values

    #calculate potential evapotranspiration for each row of table
    EPot = [E_pot(Tc[i], P[i], NETRAD[i], G[i]) for i in range(len(fluxDF))]
    EPot[EPot == 0] = np.nan 
    
    #calculate ESI
    ESI_PT = LE/EPot 

    #add as column to table
    fluxDF.insert(loc = 3, column = 'ESI PT (W/m^2)', value = ESI_PT)

    #drop any rows where ESI < 0 or greater than 1
    fluxDF = fluxDF.loc[fluxDF['ESI PT (W/m^2)'] <= upperLimit].loc[fluxDF['ESI PT (W/m^2)'] >= 0].reset_index(drop = True)
    #print(str(round(((len(fluxMonthly)-len(fluxDF))/len(fluxMonthly))*100,1)) + ' % data removed because ESI > ' + str(upperLimit) + ' or ESI < 0')
 
    return fluxDF


def minimize_b(x, y, VPD): 
    """Calculate Beta value by minimizing residuals
    Args: 
        - x (numpy array): Relative humidity from FLUXNET monthly (x-values of plot) 
        - y (numpy array): Evaporative stress index, calculated using Priestly Taylor equation and inputs from FLUXNET monthly (y-values of plot)
        - VPD (numpy array): Vapor pressure deficit in units of kPa from FLUXNET monthly
    Returns: 
        - Beta (float): Beta value, rounded to 2 decimals
    """
    #ESI parameterization
    def f(x,B):
        return x**(VPD*B)

    #Calculate residuals
    def residual(B, x, y):
        return y - f(x,B)
    
    #Minimize residuals
    def min_residual(B, x, y):
        return sum(residual(B, x, y)**2)

    #Using functions defined, minimize residuals using scipy.optimize.minimize
    B0 = 1 #initial guess at value of Beta
    method = 'L-BFGS-B' #method used to minimize a scalar
    res = optimize.minimize(min_residual, B0, method = method, args = (x, y)) 
    Beta = round(res.x[0],2)
    return Beta



def f_RH(RH,VPD,B): 
    """ESI parameterization
    Args: 
        - RH (float): Relative humidity (0-1)
        - VPD (float): Vapor pressure deficit in units of kPa
        - B (float): Beta value
    Output: 
        -ESI_B (float): Evaporative stress index, calculated using parameterization in units of W/m^2
    """
    ESI_B = RH**(VPD*B) #W/m^2
    return ESI_B 


def plotESI(RH, ESI_PT, ESI_B, Beta, title = 'Evaporative Stress Index (ESI) vs. Relative Humidity', markersize = 0.8, figPath = None): 
    """Plot actual and model ESI 
    Args: 
        - RH (np array): Relative humidity
        - ESI_PT (np array): ESI from Preistly Taylor
        - ESI_B (np array): Model ESI 
        - Beta (np array): Beta value used in model calculation 
        - title (str, optional): Title for the plot (defaults to 'Evaporative Stress Index (ESI) vs. Relative Humidity')
        - marksize (float, optional): Size of markers in scatter plot (defaults to 0.3)
        - figPath (str, optional): Path to save figure to (defaults to None-- figure not saved)
    Returns: 
        - Scatter plot displayed in notebook
    """
    #initialize figure and axes
    fig = plt.figure(figsize = (6,4))
    ax = plt.axes([0,0,1,1])

    #plot data
    ax.scatter(RH, ESI_B, s = markersize, color = 'red', label = 'f(RH) = RH^(VPD*β), β = ' + str(round(Beta,2))); #model ESI
    ax.scatter(RH, ESI_PT, s = markersize, color = 'dodgerblue', label = 'ESI tower data');

    #add descriptive info to plot
    plt.ylabel('Evaporative Stress Index')
    plt.xlabel('Relative Humidity')
    plt.xticks(np.arange(0,1+0.2, 0.2))
    plt.yticks(np.arange(0,1+0.2, 0.2))
    plt.legend(markerscale = 3, fontsize = 11.5, loc = 'upper left')
    plt.title(title, fontsize = 15)
    
    #save figure
    if figPath != None: 
        plt.savefig(figPath, bbox_inches='tight', dpi = 350)
    
    #display plot
    plt.show()


def plotResiduals(RH, residuals, title = 'Relative Humidity vs. Residuals', figPath = None): 
    """Plot residuals 
    Args: 
        - RH (np array): Relative humidity from FLUXNET monthly (x-values of plot) 
        - residuals (np array): residuals from model ESI - actual ESI 
        - title (str, optional): Title for the plot (defaults to 'Relative Humidity vs. Residuals')
        - figPath (str, optional): Path to save figure to (defaults to None-- figure not saved)
    Returns: 
        - scatter plot of residuals 
    """
    #initialize figure and axes
    fig = plt.figure(figsize = (6,4))
    ax = plt.axes([0,0,1,1])

    #plot data
    plt.scatter(RH, residuals, color = 'magenta', s = 3, label = 'Residual ')
    plt.hlines(y = 0, xmin = 0, xmax = 1, color = 'black', linestyle = '--', label = 'y = 0')

    #add descriptive info to figure
    plt.xlabel('Relative Humidity')
    plt.ylabel('Residual')
    plt.ylim(bottom = -1.05, top = 1.05)
    plt.title(title)
    plt.legend(markerscale = 3, fontsize = 11.5)

    #save figure
    if figPath != None: 
        plt.savefig(figPath, bbox_inches='tight', dpi = 350)
    
    #display figure 
    plt.show()


