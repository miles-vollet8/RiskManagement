from scholes import *
from risk import *
import plotly.graph_objects as go
import streamlit as st

st.title("Option Position Builder")
if 'positions' not in st.session_state:
    st.session_state.positions = []
if 'edit_leg_index' not in st.session_state:
    st.session_state.edit_leg_index = None

st.sidebar.header("Market Scenario")
current_price = st.sidebar.number_input("Current Underlying Price", value=100.0, step=1.0)
current_vol = st.sidebar.number_input("Current Volatility", value=0.2, step=0.01)
time_passed = st.sidebar.number_input("Time Passed (Years)", value=0.0, step=0.01)
interest_rate = st.sidebar.number_input("Interest Rate", value=RATE, step=0.001)

st.sidebar.header("Add a Leg")

leg_type = st.sidebar.radio("Leg Type", ("Option Leg", "Stock Leg"))

if leg_type == "Option Leg":
    col1, col2 = st.sidebar.columns(2)
    leg_contracts = col1.number_input("Contracts (negative for short)", value=-1, step=1)
    strike = col2.number_input("Strike Price", value=100.0, step=5.0)
    expiry = st.sidebar.number_input("Time to Expiry (Years)", value=0.5, step=0.01)
    initial_iv = st.sidebar.number_input("Initial IV", value=current_vol, step=0.01)
    entry_price = st.sidebar.number_input('Underlying Entry Price', value=current_price, step=1.00)
    option_type = st.sidebar.radio("Option Type", ("call", "put"))
elif leg_type == "Stock Leg":
    shares = st.sidebar.number_input("Shares (negative for short)", value=100, step=1)
    entry_price = st.sidebar.number_input("Entry Price", value=current_price)

if st.sidebar.button("Add Leg"):
    if leg_type == "Option Leg":
        add_leg(st.session_state.positions, leg_contracts, entry_price, strike, expiry, initial_iv, interest_rate, option_type)
        st.sidebar.success(f"Added Option Leg: {leg_contracts} contracts, Strike {strike}, Expiry {expiry}, IV {initial_iv}, Type {option_type}.")
        
    else:
        add_stock_leg(st.session_state.positions, shares, entry_price)
        st.sidebar.success(f"Added Stock Leg: {shares} shares at {entry_price}.")

# ----------------------------------------
# Display the list of legs with Remove option
# ----------------------------------------
st.subheader("Positions")
if st.session_state.positions:
    for i, leg in enumerate(st.session_state.positions):
        col1, col2, col3 = st.columns([4, 1, 1])
        with col1:
            if leg.instrument == "option":
                st.write(f"**Leg {i+1} (Option)** | Type: {leg.option_type} | Contracts: {leg.contracts} | Strike: {int(leg.strike)} | Expiry: {float(leg.time):.2f} | IV: {float(leg.initial_iv):.2f} | Entry: {leg.initial_S:.2f} | Cash Flow: {float(leg.cash_flow):.2f}")
            else:
                st.write(f"**Leg {i+1} (Stock)** | Shares: {leg.shares} | Entry Price: {float(leg.entry_price):.2f} | Cash Flow: {float(leg.cash_flow):.2f}")
        
        with col2:
            if st.button("Remove", key=f"remove_{i}"):
                st.session_state.positions.pop(i)
                st.stop()
        with col3:
            if st.button("Greeks", key=f"show_{i}"):
                if leg.instrument == "option":
                    remaining_time = max(leg.time - time_passed, 0)
                    delta, gamma, vega, theta, rho = current_leg_greeks(
                        leg, current_price, remaining_time, current_vol, interest_rate
                    )
                    st.write(f"**Leg {i+1} Greeks:**")
                    st.write(f"Delta: {delta:.2f}")
                    st.write(f"Gamma: {gamma:.2f}")
                    st.write(f"Vega: {vega:.2f}")
                    st.write(f"Theta: {theta:.2f}")
                    st.write(f"Rho: {rho:.2f}")
                else:
                    delta = 0.01*leg.shares
                    st.write(f"**Stock Leg {i+1} Greeks:**")
                    st.write(f"Delta: {delta:.2f}")
else:
    st.info("No legs added yet.")

# ----------------------------------------
# Total Position Summary: Debit/Credit and Greeks
# ----------------------------------------
has_option_legs = any(leg.instrument == "option" for leg in st.session_state.positions)

st.subheader("Total Position Summary")
if st.session_state.positions:
    total_cash = sum(leg.cash_flow for leg in st.session_state.positions)
    st.write("**Total Debit/Credit:**", f"{total_cash:.2f}")
    total_pnl = calc_pnl(st.session_state.positions, current_price, current_vol, time_passed, interest_rate)
    st.write("**Total P&L:**", f"{total_pnl:.2f}")
    if has_option_legs:
        total_delta, total_gamma, total_vega, total_theta, total_rho = position_greeks(
            st.session_state.positions, current_price, current_vol, time_passed, interest_rate
        )
        st.write("**Total Delta:**", f"{total_delta:.2f}")
        st.write("**Total Gamma:**", f"{total_gamma:.2f}")
        st.write("**Total Vega:**", f"{total_vega:.2f}")
        st.write("**Total Theta:**", f"{total_theta:.2f}")
        st.write("**Total Rho:**", f"{total_rho:.2f}")
    else:
        total_delta = sum(leg.shares * 0.01 for leg in st.session_state.positions if leg.instrument == 'stock')
        st.write("**Total Delta:**", f"{total_delta:.2f}")
        st.write("**Total Gamma:** 0.00")
        st.write("**Total Vega:** 0.00")
        st.write("**Total Theta:** 0.00")
        st.write("**Total Rho:** 0.00")
else:
    st.write("No positions to summarize.")
#st.write(f"{FLAG=}")

# ----------------------------------------
# Graph Settings
# ----------------------------------------

axis_map = {
    "Underlying": "underlying",
    "Volatility": "vol",
    "Time passed (in Years)": "time",
    "Interest Rate": "interest_rate",
    "PnL": "pnl",
    "Delta": "delta",
    "Gamma": "gamma",
    "Theta": "theta",
    "Vega": "vega",
    "Rho": "rho",
    "None": "none"
}

num_graphs = st.selectbox("Number of Graphs", [1, 2, 3], index=0)

graph_settings = []

for graph_num in range(1, num_graphs + 1):
    st.header(f"Graph Settings (Graph {graph_num})")
    col_x, col_y, col_line = st.columns(3)

    with col_x:
        x_axis = st.selectbox(
            f"X-axis (Graph {graph_num})",
            options=["Underlying", "Volatility", "Time passed (in Years)", "Interest Rate"],
            index=0,
            key=f"x_axis_{graph_num}"
        )
    with col_y:
        y_axis = st.selectbox(
            f"Y-axis (Graph {graph_num})",
            options=["PnL", "Delta", "Gamma", "Theta", "Vega", "Rho"],
            index=0,
            key=f"y_axis_{graph_num}"
        )
    with col_line:
        line_choice = st.selectbox(
            f"Extra Line (Graph {graph_num})",
            options=["None", "Underlying", "Volatility", "Time passed (in Years)", "Interest Rate"],
            index=0,
            key=f"line_choice_{graph_num}"
        )

 
    if x_axis == "Underlying":
        x_min = st.number_input(f"X-axis minimum (Graph {graph_num}, -1 for default values)", value=-1, step=1, key=f"x_min_{graph_num}")
        x_max = st.number_input(f"X-axis maximum (Graph {graph_num}, -1 for default values)", value=-1, step=1, key=f"x_max_{graph_num}")
    else:
        x_min = st.number_input(f"X-axis minimum (Graph {graph_num}, -1 for default values)", value=-1.0, step=0.01, key=f"x_min_{graph_num}")
        x_max = st.number_input(f"X-axis maximum (Graph {graph_num}, -1 for default values)", value=-1.0, step=0.01, key=f"x_max_{graph_num}")


    line_min, line_max = -1, -1
    if line_choice != "None":
        if line_choice == "Underlying":
            line_min = st.number_input(f"Line minimum (Graph {graph_num}, -1 for default values)", value=-1, step=1, key=f"line_min_{graph_num}")
            line_max = st.number_input(f"Line maximum (Graph {graph_num}, -1 for default values)", value=-1, step=1, key=f"line_max_{graph_num}")
        else:
            line_min = st.number_input(f"Line minimum (Graph {graph_num}, -1 for default values)", value=-1.0, step=0.01, key=f"line_min_{graph_num}")
            line_max = st.number_input(f"Line maximum (Graph {graph_num}, -1 for default values)", value=-1.0, step=0.01, key=f"line_max_{graph_num}")

    graph_settings.append({
        "x_axis": axis_map[x_axis],
        "y_axis": axis_map[y_axis],
        "line_choice": axis_map[line_choice],
        "x_min": x_min,
        "x_max": x_max,
        "line_values": [] if line_choice == "None" or line_min == -1 or line_max == -1 else [line_min, line_max]
    })

scenario = {
    'underlying': current_price,
    'vol': current_vol,
    'time': time_passed,
    'interest_rate': interest_rate
}

if st.button("Plot Graphs"):
    for idx, settings in enumerate(graph_settings, start=1):
        fig = plot_graph(
            st.session_state.positions,
            settings["x_axis"],
            settings["y_axis"],
            settings["line_choice"],
            scenario,
            settings["x_min"],
            settings["x_max"],
            num_points=100,
            line_values=settings["line_values"]
        )

        if fig:
            st.subheader(f"Graph {idx}")
            st.plotly_chart(fig, use_container_width=True)

