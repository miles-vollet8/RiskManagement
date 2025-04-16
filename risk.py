from scholes import *
import numpy as np
from types import SimpleNamespace
import plotly.graph_objects as go

RATE = 0.0455
def calldeltaVSunderlying(initial_S, K, time, vol, r, min_s=-1, max_s=-1,iter = 200):
    if min_s == -1:
        min_s = .7*initial_S
    if max_s == -1:
        max_s = 1.3*initial_S
    
    price_range = np.linspace(min_s, max_s, iter)
    delta_range = []
    for s in price_range:
        d1 = d1_calc(s, K, time, vol, r)
        delta = call_delta_calc(d1)
        delta_range.append(delta)

    plt.plot(price_range, delta_range)
    plt.xlabel('Underlying Price')
    plt.ylabel('Call Delta')
    plt.title('Call Delta as Underlying Changes')
    
    plt.show()

def putdeltaVSunderlying(initial_S, K, time, vol, r, min_s=-1, max_s=-1,iter = 200):
    if min_s == -1:
        min_s = .7*initial_S
    if max_s == -1:
        max_s = 1.3*initial_S
    
    price_range = np.linspace(min_s, max_s, iter)
    delta_range = []
    for s in price_range:
        d1 = d1_calc(s, K, time, vol, r)
        delta = put_delta_calc(d1)
        delta_range.append(delta)

    plt.plot(price_range, delta_range)
    plt.xlabel('Underlying Price')
    plt.ylabel('Put Delta')
    plt.title('Put Delta as Underlying Changes')
    
    plt.gca().invert_yaxis()
    plt.show()

def gammaVSunderlying(initial_S, K, time, vol, r, min_s=-1, max_s=-1,iter = 200):
    if min_s == -1:
        min_s = .7*initial_S
    if max_s == -1:
        max_s = 1.3*initial_S
    
    price_range = np.linspace(min_s, max_s, iter)
    gamma_range = []
    for s in price_range:
        d1 = d1_calc(s, K, time, vol, r)
        gamma = gamma_calc(d1, s, vol, time)
        gamma_range.append(gamma)

    plt.plot(price_range, gamma_range)
    plt.xlabel('Underlying Price')
    plt.ylabel('Gamma')
    plt.title('Gamma as Underlying Changes')
    plt.show()

def thetaVSunderlying(initial_S, K, time, vol, r, min_s=-1, max_s=-1,iter = 200):
    if min_s == -1:
        min_s = .7*initial_S
    if max_s == -1:
        max_s = 1.3*initial_S
    
    price_range = np.linspace(min_s, max_s, iter)
    theta_range = []
    for s in price_range:
        d1 = d1_calc(s, K, time, vol, r)
        d2 = d2_calc(d1, vol, time)
        theta = theta_calc(d1, d2, s, K, time, vol, r)
        theta_range.append(theta)

    plt.plot(price_range, theta_range)
    plt.xlabel('Underlying Price')
    plt.ylabel('Theta')
    plt.title('Theta as Underlying Changes')
    plt.gca().invert_yaxis()
    plt.show()

def vegaVSunderlying(initial_S, K, time, vol, r, min_s=-1, max_s=-1,iter = 200):
    if min_s == -1:
        min_s = .7*initial_S
    if max_s == -1:
        max_s = 1.3*initial_S
    
    price_range = np.linspace(min_s, max_s, iter)
    vega_range = []
    for s in price_range:
        d1 = d1_calc(s, K, time, vol, r)
        
        vega = vega_calc(d1, s, time)
        vega_range.append(vega)

    plt.plot(price_range, vega_range)
    plt.xlabel('Underlying Price')
    plt.ylabel('Vega')
    plt.title('Vega as Underlying Changes')
    plt.show()

def callrhoVStime(S, K, initial_time, vol, r, min_time=0, max_time=-1,iter = 200):
    
    if max_time == -1:
        max_time = 1.5*initial_time
        if max_time <.9:
            max_time = 1
        
    
    time_range = np.linspace(min_time, max_time, iter)
    rho_range = []
    for time in time_range:
        d1 = d1_calc(S, K, time, vol, r)
        d2 = d2_calc(d1, vol, time)
        rho = call_rho_calc(K, time, r, d2)
        rho_range.append(rho)

    plt.plot(time_range, rho_range)
    plt.xlabel('Time to Expiry')
    plt.ylabel('Call Rho')
    plt.title('Call Rho as Time Changes')
    
    plt.show()

def putrhoVStime(S, K, initial_time, vol, r, min_time=0, max_time=-1,iter = 200):
    
    if max_time == -1:
        max_time = 1.5*initial_time
        if max_time <.9:
            max_time = 1
        
    
    time_range = np.linspace(min_time, max_time, iter)
    rho_range = []
    for time in time_range:
        d1 = d1_calc(S, K, time, vol, r)
        d2 = d2_calc(d1, vol, time)
        rho = put_rho_calc(K, time, r, d2)
        rho_range.append(rho)

    plt.plot(time_range, rho_range)
    plt.xlabel('Time to Expiry')
    plt.ylabel('Put Rho')
    plt.title('Put Rho as Time Changes')
    plt.gca().invert_yaxis()
    plt.show()

#gammaVSunderlying(100, 110, .5, .2, RATE)




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
def scholes(S, K, time, vol, r, type='call'):
    
    if type == 'call':
        return call_pricer(S, K, time, vol, r)
    else: return put_pricer(S, K, time, vol, r)


def initial_leg_greeks(leg, r=RATE):
    d1 = d1_calc(leg.initial_S, leg.strike, leg.time, leg.initial_iv, r)
    d2 = d2_calc(d1, leg.initial_iv, leg.time)
    if leg.option_type == 'call':
        d1 = d1_calc(leg.inital_S, leg.strike, leg.time, leg.initial_iv, r)
        d2 = d2_calc(d1, leg.initial_iv, leg.time)
        delta = call_delta_calc(d1)
        rho = call_rho_calc(leg.strike, leg.time, r, d2)
    else: 
        delta = put_delta_calc(d1)
        rho = put_rho_calc(leg.strike, leg.time, r, d2)
    gamma = gamma_calc(d1, leg.initial_S, leg.initial_iv, leg.time)
    theta = theta_calc(d1, d2, leg.initial_S, leg.strike, leg.time, leg.initial_iv, r)
    vega = vega_calc(d1, leg.initial_S, leg.time)
    return delta, gamma, vega, theta, rho

def current_leg_greeks(leg, S, time, vol, r=RATE):
    d1 = d1_calc(S, leg.strike, time, vol, r)
    d2 = d2_calc(d1, vol, time)
    if leg.option_type == 'call':
        d1 = d1_calc(S, leg.strike, time, vol, r)
        d2 = d2_calc(d1, vol, time)
        delta = call_delta_calc(d1)
        rho = call_rho_calc(leg.strike, time, r, d2)
    else: 
        delta = put_delta_calc(d1)
        rho = put_rho_calc(leg.strike, time, r, d2)
    gamma = gamma_calc(d1, S, vol, time)
    theta = theta_calc(d1, d2, S, leg.strike, time, vol, r)
    vega = vega_calc(d1, S, time)
    if leg.contracts<0:
        delta = delta*-1
        gamma = gamma*-1
        theta = theta*-1
        vega = vega*-1
        rho = rho*-1
    return delta, gamma, vega, theta, rho

positions = []
def add_leg(positions, leg_contracts, S, K, expiry, iv, r = RATE, type='call'):
    
        
    val = scholes(S, K, expiry, iv, r, type)
    

    
    cash = -leg_contracts*100*(val)
    
    

    leg = SimpleNamespace(instrument= 'option', option_type=type, contracts=leg_contracts, strike= K, time= expiry, initial_iv=iv, initial_S=S, rate = r,cash_flow=cash)
    positions.append(leg)
    return positions

def add_stock_leg(positions, size, price):
    leg = SimpleNamespace(instrument='stock', shares=size, entry_price=price,cash_flow=-size*price)
    positions.append(leg)
    return positions
def calc_leg_pnl(leg, S=-1, vol=-1, time=0, r=-1):
    if leg.instrument == 'option':
        if S==-1:
            S = leg.initial_S
        if vol==-1:
            vol = leg.initial_iv
        if leg.time -time <0:
            time = 0
        else: time = leg.time - time
        if r == -1:
            r = leg.rate
        initial = scholes(leg.initial_S, leg.strike, leg.time, leg.initial_iv, leg.rate, leg.option_type)
        current = scholes(S, leg.strike, time, vol, r, leg.option_type)

        if leg.contracts>0:
            return leg.contracts*100*(current - initial)
        else: 
            
            print(initial)
            print(f'{current=}')
            return -leg.contracts*100*(initial - current)
    elif leg.instrument =='stock':
        if leg.shares>0:
            return leg.shares*(S -leg.entry_price)
        
        else:

            return -leg.shares*(leg.entry_price - (S))


def position_greeks(positions, S=-1, vol=-1, time=0, r=-1):
    delta = 0
    gamma = 0
    vega = 0
    theta = 0
    rho = 0

    for leg in positions:
        
        if leg.instrument == 'stock':
            
            cdelta = .01*leg.shares
            delta = delta + cdelta
            cgamma, cvega, ctheta, crho = 0, 0, 0, 0
            if leg.shares<0:
                cdelta = cdelta*-1
        elif leg.instrument =='option':
            
            current_S = leg.initial_S if S == -1 else S
            current_vol = leg.initial_iv if vol == -1 else vol
            current_r = leg.rate if r == -1 else r
            remaining_time = max(leg.time - time, 0)

            cdelta, cgamma, cvega, ctheta, crho = current_leg_greeks(
                leg, current_S, remaining_time, current_vol, current_r)
            
        delta = delta+cdelta
        gamma = gamma + cgamma
        vega = vega+cvega
        theta = theta+ ctheta
        rho = rho+crho
    

    return delta, gamma, vega, theta, rho

def calc_pnl(positions, S, vol, time, r): #leave parameters empty if they arent being affected
    pnl = 0
    for leg in positions:
        cpnl = calc_leg_pnl(leg, S, vol, time, r)
        pnl = pnl+cpnl
    return pnl

def plot_graph(positions, x_axis, y_axis, line_choice, scenario, 
               x_min, x_max, num_points=100, line_values=[]):#where scenario is the current underlying price, original testing vol, time = 0, and current interest rate
    axis_labels = {
        'underlying': 'Underlying Price',
        'vol': 'Volatility',
        'time': 'Time Passed (Years)',
        'interest_rate': 'Interest Rate'
    }
    y_labels = {
        'pnl': 'P&L',
        'delta': 'Delta',
        'gamma': 'Gamma',
        'theta': 'Theta',
        'vega': 'Vega',
        'rho': 'Rho'
    }
    defaults_x_axis = {
        'underlying': [0.5 * scenario['underlying'], 1.5 * scenario['underlying']],
        'vol': [.5*scenario['vol'], 1.5*scenario['vol']],
        'time': [0, max(1.0, scenario['time']*1.5)],#
        'interest_rate': [0.0, scenario['interest_rate'] * 2]
    }
    if x_min or x_max == -1:
        if x_min == -1:
            x_min = defaults_x_axis[x_axis][0]
        
        if x_max ==-1:
            x_max = defaults_x_axis[x_axis][1]
    if x_axis not in ['underlying', 'vol', 'time', 'interest_rate']:
        print('Invalid X axis')
        return
    if y_axis not in ['pnl', 'delta', 'gamma', 'theta', 'vega','rho']:
        print('Invalid Y axis')
        return
    if line_choice not in ['underlying', 'vol', 'time', 'interest_rate', 'none']:
        print('Invalid line choice')
        return
    if line_choice == 'x_axis':
        print('Extra line metric cannot be the same as X axis')
        return
    if not positions:
        print('No positions')
        return
    if line_choice != 'none' and not line_values:
        defaults = {
            'underlying': [0.8*scenario['underlying'], scenario['underlying'], 1.2*scenario['underlying']],
            'vol': [0.8*scenario['vol'], scenario['vol'], 1.2*scenario['vol']],
            'time': [0, 0.25, 0.5],
            'interest_rate': [0.5*scenario['interest_rate'], scenario['interest_rate'], 2*scenario['interest_rate']]
        }
        line_values = defaults[line_choice]
    elif line_choice=='none': 
        if x_axis != 'underlying':
            line_values = [scenario['underlying']]
        else: line_values = [scenario['vol']]
    
    
    x_vals = np.linspace(x_min, x_max, num_points)
    fig = go.Figure()
    i = 0
    dash_styles = ['solid', 'dash', 'dot']
    for value in line_values:
        current = scenario.copy()
        current[line_choice] = value

        y_vals = []
        for xval in x_vals:
            current[x_axis] = xval
            if y_axis == 'pnl':
                yval = calc_pnl(positions, 
                                current['underlying'], 
                                current['vol'], 
                                current['time'], 
                                current['interest_rate'])


            elif y_axis in ['delta', 'gamma', 'theta', 'vega','rho']:
                delta, gamma, vega, theta,  rho = position_greeks(positions, current['underlying'], current['vol'], current['time'], current['interest_rate'])

                if y_axis == 'delta':
                    yval = delta
                    ylabel = 'Delta'
                if y_axis == 'gamma':
                    yval = gamma
                    ylabel = 'Gamma'
                if y_axis == 'theta':
                    yval = theta
                    ylabel = 'Theta'
                if y_axis == 'vega':
                    yval = vega
                    ylabel = 'Vega'
                if y_axis == 'rho':
                    yval =rho
                    ylabel = 'Rho' 

            y_vals.append(yval)
        #linelabel = f"{line_label} at {value:.2f}"#need ifs for each possible label
        if line_choice == 'none':
            linelabel = y_labels[y_axis]
        else: linelabel = f"{axis_labels[line_choice]} = {value:.2f}"
        
        fig.add_trace(go.Scatter(x=x_vals, y=y_vals, mode='lines', name=linelabel, line=dict(dash=dash_styles[i])))
    fig.update_layout(
        title=f"{y_labels[y_axis]} vs. {axis_labels[x_axis]}",
        xaxis_title=axis_labels[x_axis],
        yaxis_title=y_labels[y_axis],
        hovermode="x unified"
    )
    return fig






# #positions = add_stock_leg(positions,50, 100)
# positions = add_leg(positions, -2, 100, 90, .1, .2, 'put')
# positions = add_leg(positions, -2, 100, 110, .1, .2)
# # delta, gamma, vega, theta, rho = position_greeks(positions, 90, [], RATE)

# # print(positions[0])
# # stock = calc_pnl(positions, 110, .2, RATE)
# # print(stock)
# # print(f'{delta=}')
# # print(f'{gamma=}')

# S_range = np.linspace(80, 120, 50)
# plot_pnl(positions, S_range)




#right below here



# S = 277
# vol = .63
# positions = []
# positions = add_leg(positions, -4, S, 330, 0.8, vol, RATE, 'call')
# positions = add_leg(positions, 4, S, 360, .8, vol, RATE, 'call')
# positions = add_leg(positions, -4, S, 240, .8, vol, RATE, 'put')
# positions = add_leg(positions, 4, S, 210, .8, vol, RATE, 'put')
# positions = add_leg(positions, -4, S, 220, .5, vol, RATE, 'put')
# positions = add_leg(positions, 4, S, 200, .5, vol, RATE, 'put')
# positions = add_leg(positions, -4, S, 320, 0.5, vol, RATE, 'call')
# positions = add_leg(positions, 4, S, 350, .5, vol, RATE, 'call')
# #positions = add_stock_leg(positions, 15, S)
# #positions = add_stock_leg(positions, 100, 100)
# scenario = {
#     'underlying': 277,
#     'vol': 0.63,
#     'time': .3,  # current time elapsed from start
#     'interest_rate': RATE
# }

# # # Interactive Plot Example
# plot_graph(
#     positions=positions,
#     x_axis='vol',
#     y_axis='pnl',
#     line_choice='time',  # You could also select 'time', 'interest_rate', or 'none'
#     scenario=scenario,x_min='none',x_max='none',
#     line_values=[.3, .5, .8]
   
# )





