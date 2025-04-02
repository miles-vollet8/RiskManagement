from math import sqrt, exp, log
import numpy as np
from scipy.stats import norm
import plotly.graph_objects as go

def d1_calc(S, K, time, vol, r): #follow call equation with input parameters
    ln = np.log(S/K) #shows how itm or otm the option is
    numerator = (ln + (r + .5*vol*vol)*time)
    denominator = (vol*np.sqrt(time))
    d1 = numerator/denominator
    return d1
def d2_calc(d1, vol, time):
    denominator = (vol*np.sqrt(time))
    d2 = d1-denominator
    return d2

def call_pricer(S, K, time, vol, r): #Where S is current price, K is strike price
    d1 = d1_calc(S, K, time, vol, r) #time is time to expiry in years, vol is volatility
    d2 = d2_calc(d1, vol, time) #and r is the risk free rate
    
    n1 = norm.cdf(d1) #delta value
    n2 = norm.cdf(d2) #probability of option expiring itm
    frac = K/(np.exp(r*time))
    call_value = S*n1 - (frac*n2)
    
    return(call_value)

def put_pricer(S, K, time, vol, r):
    call_value = call_pricer(S, K, time, vol, r) #utilize put-call parity to derive put value
    frac = K/(np.exp(r*time))
    put_value = call_value + frac - S
    return put_value

def call_delta_calc(d1):
    delta = norm.cdf(d1)
    return delta

def put_delta_calc(d1):
    delta = norm.cdf(d1) - 1
    return delta

def gamma_calc(d1, S, vol, time):
    numerator = norm.pdf(d1)
    denominator = S*vol*sqrt(time)
    gamma = numerator/denominator
    return gamma

def theta_calc(d1, d2, S, K, time, vol, r):
    term1 = (-S*norm.pdf(d1)*vol)/(2*sqrt(time))
    term2 = r*K*exp(-r*time)*norm.cdf(d2)
    theta = (term1 -term2)/365 #theta is expressed in years so divide by 365 to get daily theta
    return theta               #be aware some brokerages use 252 day convention

def vega_calc(d1, S, time):
    vega = S*norm.pdf(d1)*sqrt(time)
    return vega/100 #vega gives sensitivities to a full point(100%) of volatility change so we must divide by 100

def call_rho_calc(K, time, r, d2):
    call_rho = K * time * np.exp(-r*time) * norm.cdf(d2)
    return call_rho/100 #divide by 100 since rho is given to a full point of interest rate change

def put_rho_calc(K, time, r, d2):
    put_rho = -K * time * np.exp(-r * time)*norm.cdf(-d2)
    return put_rho/100


def call_moneyness(S,K): #call moneyness
    return S/K

def call_surface_calc(K, time, r, min_vol, max_vol, min_price, max_price): #graph with plotly
    
    colors = [(0.0, 'red'),(0.5, 'yellow'),(1.0, 'green')] #sets desired colors
    price_range = np.linspace(min_price, max_price, 50) #makes array of 50 values between min and max price
    vol_range = np.linspace(min_vol, max_vol, 50) 
    
    



    
    price_grid, vol_grid = np.meshgrid(price_range, vol_range) #maps grid
    color_vector = np.vectorize(call_moneyness)(price_grid, K) #makes vector of call moneyness for mapping colors
    call_values = np.array([call_pricer(S, K, time, vol, r) for S, vol in zip(np.ravel(price_grid), np.ravel(vol_grid))]) #makes array of call values
    call_values = call_values.reshape(price_grid.shape)
    
    call_surface = go.Figure(data=[go.Surface(z=call_values, x=price_grid, y=vol_grid, surfacecolor = color_vector, #create call surface
                                              colorscale=colors, cmin = 0.7, cmax = 1.3, cmid=1,
                                              colorbar=dict(
                                                tickvals=[0.7, 1.0, 1.3],  # Correspond to OTM, ATM, ITM
                                                ticktext=["OTM", "ATM", "ITM"],  # Labels for colors
                                                          
                                            ),
                                              hovertemplate=(
            "Asset Price: %{x}<br>" + #improves readability of graphs
            "Volatility: %{y}<br>" +
            "Call Value: %{z}<extra></extra>"))])
    call_surface.update_layout(
        
        scene=dict(
            xaxis_title="Underlying Price (S)", #axis labels
            yaxis_title="Volatility (σ)",
            zaxis_title="Call Option Value"
        )
    )
  
    return call_surface

def put_moneyness(S,K): #put moneyness is call reciprocal
    return K/S

def put_surface_calc(K, time, r, min_vol, max_vol, min_price, max_price):
    colors = [(0.0, 'red'),(0.5, 'yellow'),(1.0, 'green')]

    price_range = np.linspace(min_price, max_price, 50)
    vol_range = np.linspace(min_vol, max_vol, 50)
    
    
    price_grid, vol_grid = np.meshgrid(price_range, vol_range)
    color_vector = np.vectorize(put_moneyness)(price_grid, K)
    put_values = np.array([put_pricer(S, K, time, vol, r) for S, vol in zip(np.ravel(price_grid), np.ravel(vol_grid))])#prices puts instead of calls
    put_values = put_values.reshape(price_grid.shape)
    put_surface = go.Figure(data=[go.Surface(z=put_values, x=price_grid, y=vol_grid, surfacecolor = color_vector, 
                                             colorscale=colors, cmin = 0.7, cmax = 1.3, cmid=1,
                                             colorbar=dict(
                                                tickvals=[0.7, 1.0, 1.3],
                                                ticktext=["OTM", "ATM", "ITM"],
                                            ),
                                             hovertemplate=(
            "Asset Price: %{x}<br>" +
            "Volatility: %{y}<br>" +
            "Put Value: %{z}<extra></extra>"))])
    put_surface.update_layout(
        
        scene=dict(
            xaxis_title="Underlying Price (S)",
            yaxis_title="Volatility (σ)",
            zaxis_title="Put Option Value",
            
        )
    )
    
    return put_surface

