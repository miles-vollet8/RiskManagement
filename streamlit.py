import streamlit as st
from scholes import *
st.title("Black-Scholes Model")
st.write("Find european option prices and sensitivities with the Black-Scholes formula")

st.markdown("""
<style>
.metric-container {
    display: flex;
    justify-content: center; /* centers text */
    padding: 8px;  /* gives some space around text */ 
    width: auto;     margin: 0 auto; /* center container */
    }
.metric-call{
            background-color: #90EE90; 
            color: black;
            margin-right: 10px;
            border-radius: 10px;
            }    
.metric-put{
            background-color: #FF6347;
            color: black;
            border-radius: 10px;
            }

.metric-label{
            font-size: 1rem;

            } 
.metric-value{
            font-size: 1.5rem; /* makes option values more prominent */
            font-weight: bold; 
            }                                   
            </style>
            """, unsafe_allow_html=True)





with st.sidebar:
    st.title("Black-Scholes Model")
    st.markdown(
    '''
    Created By:
    <a href="https://www.linkedin.com/in/miles-vollet-9044a1314/" target="_blank"> 
        <img src="https://cdn-icons-png.flaticon.com/512/174/174857.png" alt="LinkedIn" width="20" style="vertical-align: middle; margin-right: 5px;">
    </a>
    <a href="https://www.linkedin.com/in/miles-vollet-9044a1314/" target="_blank">
        Miles Vollet
    </a> 
    ''', #linkedin plug
    unsafe_allow_html=True)

    S = st.sidebar.number_input("Current Asset Price(S)", min_value=0.0, value=100.00, step=1.0) #allows user inputs
    K = st.sidebar.number_input("Strike Price(K)", min_value=0.0, value=100.00, step=1.0)
    time = st.sidebar.number_input("Time to Expiry (in years)", min_value=0.0, value=1.00, step=0.01)
    vol = st.sidebar.number_input("Volatility(σ)", min_value=0.0, value=.25, step=0.01)
    r= st.sidebar.number_input("Risk-Free Interest Rate", value=.05, step=0.01)

    st.markdown("---")
    surface_gen = st.button("Generate Surface")
    
    min_price = st.number_input('Minimum Asset Price', min_value=0.01, value=S*.75, step=1.0)
    max_price = st.number_input("Maximum Asset Price", min_value=0.01, value=S*1.25, step=1.0)
    min_vol = st.number_input("Minimum Volatility", min_value=0.0, value=vol*.25, step=0.01)
    max_vol = st.number_input("Maximum Volatility", min_value=0.01, value=vol*1.75, step = 0.01)
    price_range = np.linspace(min_price, max_price, 100)
    vol_range = np.linspace(min_vol, max_vol, 100)


call_val = call_pricer(S,K,time,vol,r) #calculations
d1 = d1_calc(S, K, time, vol, r)
d2 = d2_calc(d1, vol, time)
call_delta = call_delta_calc(d1)
put_delta = put_delta_calc(d1)
gamma = gamma_calc(d1, S, vol, time)
theta = theta_calc(d1, d2, S, K, time, vol, r)
vega = vega_calc(d1, S, time)
call_rho = call_rho_calc(K, time, r, d2)
put_rho = put_rho_calc(K, time, r, d2)
put_val = put_pricer(S,K,time,vol,r)
col1, col2 = st.columns([1,1]) #sets up call and put columns

with col1:
    st.markdown("""
    <style>
    /* Adjust column widths */
    .css-1lcbmhc { width: 48% !important; }
    .css-1q9ozre { width: 48% !important; }
    .css-1kyxreq { margin-left: 1% !important; }
    </style>
    """, unsafe_allow_html=True)
    st.markdown(f"""
                <div class="metric-container metric-call">
    <div>
        <div class="metric-label">Call Value</div>
        <div class="metric-value">${call_val:.2f}</div>
    </div>
    </div>
    """, unsafe_allow_html=True)
    st.write(f"**Delta(Δ)**: {call_delta:.4f}") #displays greeks
    st.write(f"**Gamma(Γ)**: {gamma:.4f}")
    st.write(f"**Theta(Θ)**: {theta:.4f}")
    st.write(f"**Vega(ν)**: {vega:.4f}")
    st.write(f"**Rho(ρ)**: {call_rho:.4f}")
    
    
with col2:
    st.markdown("""
    <style>
    /* Adjust column widths */
    .css-1lcbmhc { width: 48% !important; }
    .css-1q9ozre { width: 48% !important; }
    .css-1kyxreq { margin-left: 1% !important; }
    </style>
    """, unsafe_allow_html=True)
    st.markdown(f"""
        <div class="metric-container metric-put">
            <div>
                <div class="metric-label">Put Value</div>
                <div class="metric-value">${put_val:.2f}</div>
            </div>
        </div>
    """, unsafe_allow_html=True)
    st.write(f"**Delta(Δ)**: {put_delta:.4f}")
    st.write(f"**Gamma(Γ)**: {gamma:.4f}")
    st.write(f"**Theta(Θ)**: {theta:.4f}")
    st.write(f"**Vega(ν)**: {vega:.4f}")
    st.write(f"**Rho(ρ)**: {put_rho:.4f}")
    
            
st.markdown("---")
st.title("Option Price Interactive Surface")
st.write("See how Option Prices(Z) fluctuate depending on Volatility(Y) and Underlying Asset changes(X)")
col1, col2 = st.columns([1,1])
with col1:
    st.subheader("Call Option")
    call_surface = call_surface_calc(K, time, r, min_vol, max_vol, min_price, max_price)
    st.plotly_chart(call_surface) #plot call surface

with col2:
    st.subheader("Put Option")
    put_surface = put_surface_calc(K, time, r, min_vol, max_vol, min_price, max_price)
    st.plotly_chart(put_surface) #plot put surface





