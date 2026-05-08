# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo",
#     "numpy==2.4.4",
#     "altair==6.1.0",
#     "pandas==3.0.2",
#     "pyarrow==24.0.0",
# ]
# ///

import marimo

__generated_with = "0.23.5"
app = marimo.App(
    width="medium",
    css_file="/usr/local/_marimo/custom.css",
    auto_download=["html"],
)


@app.cell(hide_code=True)
def imports():
    # Bibliotheken
    import marimo as mo
    import numpy as np
    import altair as alt
    import pandas as pd

    return alt, mo, np, pd


@app.cell(hide_code=True)
def titel(mo):
    # Titelkaart
    mo.md("""
    # Methylfenidaat — Plasmaverloop

    Dit notebook modelleert hoe methylfenidaat (MPH) zich gedurende de dag in je bloed
    verspreidt, op basis van jouw lichaamskenmerken. Je kiest een doseerschema, stelt je
    ontwaaktijd in, en de grafiek toont direct de verwachte plasmaconcentratie, piekmomenten
    en het hypothetische rebound-moment.

    **Functies:** gepersonaliseerde PK-curves op basis van gewicht, lengte en metabolisme —
    vier doseerschema’s met verstelbare timers — rebound-detectie in de grafiek —
    woordelijke samenvatting van je dagverloop.
    """)
    return


@app.cell(hide_code=True)
def lichaamskenmerken_defaults(mo):
    # Lichaamskenmerken — invoerwidgets
    sl_weight = mo.ui.slider(start=40, stop=140, step=1, value=75,
        label="Lichaamsgewicht (kg)", show_value=True)
    sl_height = mo.ui.slider(start=140, stop=210, step=1, value=174,
        label="Lichaamslengte (cm)", show_value=True)
    sex_radio = mo.ui.radio(options={"Man": "M", "Vrouw": "V"},
        value="Man", label="Geslacht", inline=True)
    sl_metab = mo.ui.slider(steps=[-0.5, -0.3, 0.0, 0.3, 0.5], value=0.0,
        label="Metabolisme (CES1): langzaam ← gemiddeld → snel", show_value=True)
    return sex_radio, sl_height, sl_metab, sl_weight


@app.cell(hide_code=True)
def pk_parameters(np, sex_radio, sl_height, sl_metab, sl_weight):
    # PK-parameters — berekend uit lichaamskenmerken
    _W   = float(sl_weight.value)
    _H   = float(sl_height.value)
    _sex = sex_radio.value
    _met = float(sl_metab.value)
    # Body Mass Index — nodig voor LBM-formule
    _BMI = _W / (_H / 100) ** 2
    # Lean Body Mass via Janmahasatian (2005) — basis voor distributievolume Vd
    _LBM = (9270*_W)/(6680+216*_BMI) if _sex=="M" else (9270*_W)/(8780+244*_BMI)
    _LBM_ref   = 56.0
    # Allometrische klaring-schaling: CL ∝ gewicht^0.75 (Lyauk et al. 2016)
    _CL_factor = (_W / 70.0) ** 0.75
    # CES1-metabolisme: exponent op slider zodat ±0.5 ≈ factor 1.65×
    _mf        = float(np.exp(_met))
    # Vd schaalt met LBM t.o.v. referentie (56 kg, gemiddelde man)
    _Vd_factor = _LBM / _LBM_ref
    # Netto eliminatie-schaling: ke_scale > 1 → kortere halfwaardetijd
    _ke_scale  = (_CL_factor * _mf) / _Vd_factor
    pk_params = {
        "ke_scale":     float(_ke_scale),
        "scale_factor": float(1.0 / _Vd_factor),
        "CL_factor":    float(_CL_factor),
        "metab_factor": float(_mf),
        "LBM":          float(_LBM),
        "BMI":          float(_BMI),
        # t½ = ln(2) / ke; hogere ke_scale → kortere halfwaardetijd
        "t_half_ir":    round(float(np.log(2)/(0.277*_ke_scale)), 1),
        "t_half_er":    round(float(np.log(2)/(0.198*_ke_scale)), 1),
    }
    return (pk_params,)


@app.cell(hide_code=True)
def pk_functies(np, pk_params):
    # PK-functies — IR- en ER-concentratiecurven
    # Basis PK-constanten (populatie-gemiddelden, fixed)
    IR_KA=3.5; IR_KE0=0.277; ER_KA=0.38; ER_KE0=0.198; SCALE0=0.40
    # Schaal eliminatie (ke) en concentratie (SCALE) met lichaamsparameters
    _ks=pk_params["ke_scale"]; _sf=pk_params["scale_factor"]
    IR_KE=IR_KE0*_ks; ER_KE=ER_KE0*_ks; SCALE=SCALE0*_sf

    # AUC-genormaliseerde Bateman-functie: integraal = 1 voor elk ke/ka-paar
    def bateman_normed(dt, ka, ke):
        raw=np.exp(-ke*dt)-np.exp(-ka*dt)
        return np.maximum(raw/(1.0/ke-1.0/ka), 0.0)

    # IR-curve: enkelvoudige Bateman-piek, actief gedurende 14 h na inname
    def ir_curve(t, dose, t0):
        dt=t-t0
        return np.where((dt>=0)&(dt<=14), dose*SCALE*bateman_normed(np.maximum(dt,0),IR_KA,IR_KE), 0.0)

    # ER-curve: 22% snelle IR-coating + 78% trage osmotische pomprelease
    def er_curve(t, dose, t0):
        dt=t-t0
        c=dose*0.22*SCALE*bateman_normed(np.maximum(dt,0),IR_KA,IR_KE) + \
          dose*0.78*SCALE*bateman_normed(np.maximum(dt,0),ER_KA,ER_KE)
        return np.where((dt>=0)&(dt<=16), c, 0.0)

    return er_curve, ir_curve


@app.cell(hide_code=True)
def y_schaal(np):
    # Grafiek — Y-as referentieschaal (referentielichaam 70 kg)
    # Y_MAX wordt berekend met vaste referentie-PK-constanten — onafhankelijk van
    # lichaamskenmerken. Zo worden wijzigingen in gewicht/metabolisme zichtbaar als
    # hoogte-verandering in de grafiek (hogere/lagere curves t.o.v. de vaste as).
    _KA_IR, _KE_IR = 3.5, 0.277
    _KA_ER, _KE_ER = 0.38, 0.198
    _SCALE = 0.40

    def _bat(dt, ka, ke):
        raw = np.exp(-ke * np.maximum(dt, 0)) - np.exp(-ka * np.maximum(dt, 0))
        return np.maximum(raw / (1/ke - 1/ka), 0.0)

    def _ir(t, dose, t0):
        dt = t - t0
        return np.where((dt >= 0) & (dt <= 14), dose * _SCALE * _bat(dt, _KA_IR, _KE_IR), 0.0)

    def _er(t, dose, t0):
        dt = t - t0
        c = (dose * 0.22 * _SCALE * _bat(dt, _KA_IR, _KE_IR) +
             dose * 0.78 * _SCALE * _bat(dt, _KA_ER, _KE_ER))
        return np.where((dt >= 0) & (dt <= 16), c, 0.0)

    _t = np.linspace(6.0, 23.0, 800)
    def _som(specs): return sum(fn(_t, mg, 7.0 + dh) for fn, mg, dh in specs)
    # Maximum over alle vier schema’s, plus 18% marge
    Y_MAX = float(max(
        _som([(_ir, 10, 0), (_ir, 10, 4), (_ir, 5, 8)]).max(),
        _som([(_er, 27, 0)]).max(),
        _som([(_er, 27, 0), (_ir, 10, 6)]).max(),
        _som([(_ir, 10, 0), (_er, 27, 0.5)]).max(),
    )) * 1.18
    return (Y_MAX,)


@app.cell(hide_code=True)
def schema_sliders(mo):
    # Schema’s — tijdoffset-sliders
    def _s(label, default, min_h=0.5):
        return mo.ui.slider(steps=[h/2 for h in range(int(min_h*2),29)],
                            value=default, label=label, show_value=True)
    s1_d2=_s("2e dosis 10 mg IR — uur na ontwaaktijd (standaard 4 h)", 4.0)
    s1_d3=_s("3e dosis 5 mg IR — uur na ontwaaktijd (standaard 8 h)",  8.0)
    s3_ir=_s("10 mg IR booster — uur na ontwaaktijd (standaard 6 h)",  6.0)
    s4_er=_s("27 mg Concerta ER — uur na ontwaaktijd (standaard 0.5 h)", 0.5)
    return s1_d2, s1_d3, s3_ir, s4_er


@app.cell(hide_code=True)
def ontwaaktijd_invoer(mo):
    # Ontwaaktijd — invoer
    wake_input = mo.ui.text(value="07:00", label="Ontwaaktijd / 1e dosis (HH:MM)",
                            placeholder="07:00")
    return (wake_input,)


@app.cell(hide_code=True)
def sectie_grafiek(mo):
    # Sectie — Grafiek
    mo.md("""
    ## Grafiek
    """)
    return


@app.cell(hide_code=True)
def grafiek(
    Y_MAX,
    alt,
    er_curve,
    ir_curve,
    mo,
    np,
    pd,
    s1_d2,
    s1_d3,
    s3_ir,
    s4_er,
    schema_tabs_ui,
    wake_input,
):
    # Grafiek — concentratieverloop over de dag

    # Kleurenpalet
    _C_IR10 = "#1565C0"   # donkerblauw — IR 10 mg
    _C_IR5  = "#42A5F5"   # lichtblauw  — IR 5 mg
    _C_ER27 = "#F57F17"   # geel/amber  — Concerta ER 27 mg
    _C_TOT  = "#2E7D32"   # donkergroen — totale curve
    _C_REB  = "#C62828"   # donkerrood  — rebound

    # Ontwaaktijd parsen naar uren (bijv. "07:30" → 7.5)
    try:
        _h,_m=wake_input.value.strip().split(":")
        wake=int(_h)+int(_m)/60
    except Exception:
        wake=7.0

    # Actief schema bepalen en dosis-lijst opbouwen
    active=schema_tabs_ui.value

    if active=="Schema 1 — Alleen ER":
        doses=[{"mg":27,"fn":er_curve,"t0":wake,"label":"27 mg Concerta ER","color":_C_ER27}]
    elif active=="Schema 2 — Alleen IR":
        doses=[
            {"mg":10,"fn":ir_curve,"t0":wake,            "label":"10 mg IR",        "color":_C_IR10},
            {"mg":10,"fn":ir_curve,"t0":wake+s1_d2.value, "label":"10 mg IR",        "color":_C_IR10},
            {"mg": 5,"fn":ir_curve,"t0":wake+s1_d3.value, "label":"5 mg IR (half)",  "color":_C_IR5},
        ]
    elif active=="Schema 3 — ER + IR":
        doses=[
            {"mg":27,"fn":er_curve,"t0":wake,            "label":"27 mg Concerta ER","color":_C_ER27},
            {"mg":10,"fn":ir_curve,"t0":wake+s3_ir.value, "label":"10 mg IR",         "color":_C_IR10},
        ]
    else:
        doses=[
            {"mg":10,"fn":ir_curve,"t0":wake,            "label":"10 mg IR",          "color":_C_IR10},
            {"mg":27,"fn":er_curve,"t0":wake+s4_er.value, "label":"27 mg Concerta ER","color":_C_ER27},
        ]

    # Tijdas: van 24 min vóór ontwaak tot 15:36 erna, 700 punten
    t_arr=np.linspace(wake-0.4, wake+15.6, 700)

    # Hulpfunctie: uren naar HH:MM string
    def fmt(h):
        H=int(h)%24; M=round((h%1)*60)
        if M==60: H,M=(H+1)%24,0
        return f"{H:02d}:{M:02d}"

    # X-as ticks: elk heel uur, omgezet naar HH:MM via Vega-expressie
    tick_vals=list(range(int(wake),int(wake+16)))
    lmap_js="{"+", ".join(f"{v}: '{fmt(v)}'" for v in tick_vals)+"}"

    # Bereken individuele curves en sommeer tot totale curve
    rows,total=[],np.zeros_like(t_arr)
    for d in doses:
        c=d["fn"](t_arr,d["mg"],d["t0"]); total+=c
        for ti,ci in zip(t_arr,c):
            rows.append({"t":ti,"Tijd":fmt(ti),"ng_mL":round(float(ci),4),"Dosis":d["label"]})

    df=pd.DataFrame(rows)
    # Dosis-label "Totaal" zodat de kleurschaal de totale lijn meeneemt in de legenda
    df_total=pd.DataFrame([
        {"t":float(ti),"Tijd":fmt(float(ti)),"ng_mL":round(float(ci),4),"Dosis":"Totaal"}
        for ti,ci in zip(t_arr,total)
    ])

    # Rebound-detectie: eerste punt na dagpiek waar totale curve < 25% van piek
    peak_val  = float(total.max())
    threshold = peak_val*0.25
    tmax_idx  = int(np.argmax(total))
    tmax_t    = t_arr[tmax_idx]
    rb_mask   = (t_arr>tmax_t)&(total<threshold)
    rebound_t = float(t_arr[np.argmax(rb_mask)]) if rb_mask.any() else None

    # Zone achter reboundmoment (rode achtergrondvlak)
    df_rb_zone=pd.DataFrame([{"t":float(ti),"ng_mL":round(float(ci),4)}
                              for ti,ci in zip(t_arr,total)
                              if rebound_t is not None and ti>=rebound_t])
    # Horizontale drempellijn van dagpiek tot einde tijdvenster
    df_thr=pd.DataFrame([{"t":float(tmax_t),"ng_mL":threshold},
                          {"t":float(t_arr[-1]),"ng_mL":threshold}])

    # Piekmomenten per afzonderlijke dosis (voor woordelijke samenvatting)
    dose_peaks=[]
    for d in doses:
        _c=d["fn"](t_arr,d["mg"],d["t0"])
        if _c.max()>0.01:
            dose_peaks.append({"label":d["label"],"mg":d["mg"],"t0":d["t0"],
                                "tpeak":float(t_arr[np.argmax(_c)]),"peak":float(_c.max())})

    # Kleurschaal: individuele doses + Totaal + Rebound — alle vijf zichtbaar in legenda
    labels=list(dict.fromkeys(d["label"] for d in doses))
    colors=[next(d["color"] for d in doses if d["label"]==l) for l in labels]
    labels_full = labels + ["Totaal", "Rebound"]
    colors_full = colors + [_C_TOT, _C_REB]
    _cscale = alt.Scale(domain=labels_full, range=colors_full)

    # Y-domein: vaste referentie als minimum; groeit mee als huidige params hoger uitkomen
    _y_domain = max(Y_MAX, float(total.max()) * 1.05)

    # Altair-encodings
    x=alt.X("t:Q",scale=alt.Scale(domain=[wake-0.4,wake+15.6]),
            axis=alt.Axis(values=tick_vals,labelExpr=f"({lmap_js})[datum.value] || ''",
                          title=None,grid=True,gridColor="#e8e8e8",labelColor="#555",tickColor="#bbb"))
    y=alt.Y("ng_mL:Q",scale=alt.Scale(domain=[0,_y_domain]),
            axis=alt.Axis(title="ng/mL (gepersonaliseerd)",grid=True,gridColor="#e8e8e8",labelColor="#555"))
    y_tot=alt.Y("ng_mL:Q",scale=alt.Scale(domain=[0,_y_domain]))
    # Hoofdlegenda met alle vijf series (Altair toont ook domain-waarden die niet in df zitten)
    col=alt.Color("Dosis:N",scale=_cscale,legend=alt.Legend(title=None,orient="top-left"))

    # Gekleurde vlakken per dosis (lichte vulling)
    areas=alt.Chart(df).mark_area(opacity=0.13,interpolate="monotone").encode(x=x,y=y,color=col)
    # Gekleurde lijnen per dosis (geen dubbele legenda)
    lines=alt.Chart(df).mark_line(strokeWidth=2.2,interpolate="monotone").encode(
        x=x,y=y,color=alt.Color("Dosis:N",scale=_cscale,legend=None))
    # Totale curve: groene stippellijn, kleur uit gedeelde schaal
    total_line=alt.Chart(df_total).mark_line(strokeWidth=2.8,strokeDash=[7,3],
        interpolate="monotone").encode(
        x=x,y=y_tot,color=alt.Color("Dosis:N",scale=_cscale,legend=None))

    # Rode zone achter reboundmoment
    rb_area=(alt.Chart(df_rb_zone).mark_area(opacity=0.18,color=_C_REB,interpolate="monotone")
             .encode(x=x,y=y_tot)) if len(df_rb_zone)>0 \
        else alt.Chart(pd.DataFrame({"t":[],"ng_mL":[]})).mark_point()

    # Horizontale drempellijn (25%-niveau)
    thr_line=alt.Chart(df_thr).mark_line(strokeDash=[4,4],strokeWidth=1.5,
        color=_C_REB,opacity=0.7).encode(
        x=alt.X("t:Q"),y=alt.Y("ng_mL:Q",scale=alt.Scale(domain=[0,_y_domain])))

    # Verticale stippellijnen op inname-tijdstippen
    vlines_df=pd.DataFrame([{"t":float(d["t0"]),"label":f"▼ {d['mg']} mg  ({fmt(d['t0'])})","color":d["color"]}
                             for d in doses])
    dose_rules=alt.Chart(vlines_df).mark_rule(strokeDash=[3,3],strokeWidth=1.5,opacity=0.45).encode(
        x=alt.X("t:Q"),color=alt.Color("color:N",scale=None),tooltip=alt.Tooltip("label:N"))

    # Verticale rebound-lijn met kleur uit gedeelde schaal
    rb_layers=[]
    if rebound_t is not None:
        _df_rb_rule=pd.DataFrame([{"t":rebound_t,"label":f"Rebound ~{fmt(rebound_t)}","Dosis":"Rebound"}])
        rb_layers.append(
            alt.Chart(_df_rb_rule)
            .mark_rule(strokeWidth=2,opacity=0.8)
            .encode(
                x=alt.X("t:Q"),
                color=alt.Color("Dosis:N",scale=_cscale,legend=None),
                tooltip=alt.Tooltip("label:N")))

    # Interactieve cursor + dot bij hover over totale curve
    nearest=alt.selection_point(nearest=True,on="mouseover",fields=["t"],empty=False)
    selectors=(alt.Chart(df_total).mark_point(opacity=0)
        .encode(x=x,tooltip=[alt.Tooltip("Tijd:N",title="Tijd"),
                              alt.Tooltip("ng_mL:Q",title="Totaal (ng/mL)",format=".2f")])
        .add_params(nearest))
    cursor=(alt.Chart(df_total).mark_rule(color="#aaa",strokeWidth=1)
        .encode(x=x).transform_filter(nearest))
    dot=(alt.Chart(df_total).mark_point(filled=True,size=70,color=_C_TOT)
        .encode(x=x,y=y_tot,opacity=alt.condition(nearest,alt.value(1),alt.value(0)),
                tooltip=[alt.Tooltip("Tijd:N",title="Tijd"),
                         alt.Tooltip("ng_mL:Q",title="Totaal ng/mL",format=".2f")]))

    chart=alt.layer(rb_area,thr_line,areas,lines,total_line,dose_rules,selectors,cursor,dot,
                    *rb_layers
    ).properties(width="container",height=340,title=active).configure_view(
        strokeOpacity=0,fill="white").configure(background="white",font="system-ui, sans-serif")

    mo.ui.altair_chart(chart)
    return dose_peaks, fmt, peak_val, rebound_t, threshold


@app.cell(hide_code=True)
def grafiek_beschrijving(
    dose_peaks,
    fmt,
    mo,
    peak_val,
    pk_params,
    rebound_t,
    schema_tabs_ui,
    threshold,
):
    # Grafiek — woordelijke samenvatting
    _active = schema_tabs_ui.value
    _tir    = pk_params["t_half_ir"]
    _ter    = pk_params["t_half_er"]

    # Inname-tijden
    _innamen = [f"{fmt(d['t0'])} ({d['label']}, {d['mg']} mg)" for d in dose_peaks]
    _innamen_str = ", ".join(_innamen) if _innamen else "—"

    # Piekmomenten
    _pieken = [f"{fmt(d['tpeak'])} ({d['label']}: {d['peak']:.2f} ng/mL)" for d in dose_peaks]
    _pieken_str = ", ".join(_pieken) if _pieken else "—"

    _rb_str = (f"~{fmt(rebound_t)} (concentratie daalt onder {threshold:.2f} ng/mL = 25% van piek)"
               if rebound_t is not None else "niet verwacht binnen het dagvenster")

    mo.accordion({
        "Woordelijke samenvatting van dit schema (klik om te lezen)": mo.md(f"""
    **Schema:** {_active}

    **Innamen:** {_innamen_str}

    **Piekmomenten per dosis:** {_pieken_str}

    **Dagpiek (totale curve):** {peak_val:.2f} ng/mL

    **Werkzame drempel (25% van piek):** {threshold:.2f} ng/mL

    **Verwacht rebound-moment:** {_rb_str}

    **Halfwaardetijden:** IR = {_tir} h, Concerta ER = {_ter} h

    De groene stippellijn toont de totale plasmaspiegel gedurende de dag.
    De gekleurde vlakken tonen de bijdrage van elke afzonderlijke dosis.
    De rode zone markeert de periode vanaf het hypothetische rebound-moment,
    wanneer de concentratie onder de werkzame drempel daalt.
        """)
    })
    return


@app.cell(hide_code=True)
def sectie_schema(mo):
    # Sectie — Schema & tijdstip
    mo.md("""
    ## Schema & tijdstip
    """)
    return


@app.cell
def schema_invoer(wake_input):
    # Ontwaaktijd — weergave
    wake_input
    return


@app.cell(hide_code=True)
def schema_tabs_cel(mo, s1_d2, s1_d3, s3_ir, s4_er, wake_input):
    # Schema-tabs — selectie en wekker-tijden

    # Helperfuncties voor tijdnotatie
    def _parse(v):
        try:
            h,m=v.strip().split(":"); return int(h)+int(m)/60
        except Exception: return 7.0

    def _hhmm(h):
        H=int(h)%24; M=int(round((h%1)*60))
        if M==60: H,M=(H+1)%24,0
        return f"{H:02d}:{M:02d}"

    _wake=_parse(wake_input.value)

    # Wekker-callout: geeft absolute kloktijden per dosis
    def _alarm(doses):
        lines=["**Wekker zetten op:**"]
        for lbl,dh in doses: lines.append(f"- **{_hhmm(_wake+dh)}** — {lbl}")
        return mo.callout(mo.md("\n".join(lines)),kind="success")

    schema_tabs_ui=mo.ui.tabs({
        "Schema 1 — Alleen ER": mo.vstack([
            mo.md("**27 mg Concerta ER** — referentieschema, een pil per dag."),
            _alarm([("27 mg Concerta ER",0)]),
        ],gap="0.5rem"),
        "Schema 2 — Alleen IR": mo.vstack([
            mo.md("**10 mg IR + 10 mg IR + 5 mg IR** — volledige controle, geen ER."),
            mo.md(f"1e dosis bij ontwaaktijd ({_hhmm(_wake)})"),
            s1_d2, s1_d3,
            _alarm([("10 mg IR (1e dosis)",0),("10 mg IR (2e dosis)",s1_d2.value),("5 mg IR (3e dosis)",s1_d3.value)]),
        ],gap="0.5rem"),
        "Schema 3 — ER + IR": mo.vstack([
            mo.md("**27 mg Concerta ER** bij ontwaaktijd + **10 mg IR** booster later."),
            s3_ir,
            _alarm([("27 mg Concerta ER",0),("10 mg IR booster",s3_ir.value)]),
        ],gap="0.5rem"),
        "Schema 4 — IR + ER": mo.vstack([
            mo.md("**10 mg IR** bij ontwaaktijd + **27 mg Concerta ER** iets later."),
            s4_er,
            _alarm([("10 mg IR (snelle start)",0),("27 mg Concerta ER",s4_er.value)]),
        ],gap="0.5rem"),
    })
    schema_tabs_ui
    return (schema_tabs_ui,)


@app.cell(hide_code=True)
def sectie_lichaam(mo):
    # Sectie — Lichaamskenmerken
    mo.md("""
    ## Lichaamskenmerken
    """)
    return


@app.cell(hide_code=True)
def lichaamskenmerken_ui(
    mo,
    pk_params,
    sex_radio,
    sl_height,
    sl_metab,
    sl_weight,
):
    # Lichaamskenmerken — weergave en PK-tabel
    mo.vstack([
        mo.hstack([sl_weight, sl_height], gap="2rem", justify="start"),
        sex_radio,
        sl_metab,
        mo.accordion({
            "Wat doen deze parameters? (klik om uit te klappen)": mo.lazy(lambda: mo.md("""
    **Gewicht → klaring (CL)**
    Zwaardere mensen clearen MPH sneller via CES1 in de lever.
    Klaring schaalt met gewicht^0.75 (allometrische wet, Lyauk et al. 2016).
    Hogere CL = lagere piek, kortere duur.

    **Lengte + gewicht → lean body mass (LBM)**
    MPH lost slecht op in vet (laag log P). Het distributievolume (Vd) hangt
    af van vetvrije massa, berekend via de Janmahasatian-formule (2005).

    **Metabolisme-slider (CES1-activiteit)**
    CES1-enzymactiviteit varieert ~3-4x door genetica — niet meetbaar zonder lab.
    Stel dit in als je merkt dat de medicatie sneller of langzamer uitwerkt
    dan verwacht. Links = langzame metaboliseerder (PM), rechts = snelle (UM).
            """))
        }),
        mo.accordion({
            f"Berekende PK-parameters — t½ IR: {pk_params['t_half_ir']} h, t½ ER: {pk_params['t_half_er']} h":
            mo.lazy(lambda: mo.md(f"""
    | Parameter | Waarde | Effect |
    |-----------|--------|--------|
    | BMI | {pk_params['BMI']:.1f} kg/m² | — |
    | Lean Body Mass | {pk_params['LBM']:.1f} kg | Basis voor Vd |
    | Klaring-factor | {pk_params['CL_factor']:.2f}× | {'hoger gewicht → snellere afbraak' if pk_params['CL_factor']>1 else 'lager gewicht → langzamere afbraak'} |
    | Metabolisme-factor | {pk_params['metab_factor']:.2f}× | {'sneller CES1' if pk_params['metab_factor']>1 else 'langzamer CES1'} |
    | Netto ke-schaling | {pk_params['ke_scale']:.2f}× | t½ IR: {pk_params['t_half_ir']} h, t½ ER: {pk_params['t_half_er']} h |
    | Concentratie-schaling | {pk_params['scale_factor']:.2f}× | {'lagere piek' if pk_params['scale_factor']<1 else 'hogere piek'} |
            """))
        }),
    ], gap="0.6rem")
    return


@app.cell(hide_code=True)
def sectie_rebound(mo):
    # Sectie — Rebound
    mo.md("""
    ## Rebound
    """)
    return


@app.cell(hide_code=True)
def rebound_info(fmt, mo, peak_val, rebound_t, threshold):
    # Rebound — callout en uitleg
    _callout = (
        mo.callout(mo.md(
            f"**Hypothetisch rebound-moment: ~{fmt(rebound_t)}**  \n"
            f"Concentratie daalt onder de werkzame drempel "
            f"({threshold:.2f} ng/mL = 25% van dagpiek {peak_val:.2f} ng/mL).  \n"
            f"Rode zone in de grafiek = periode met verhoogde kans op reboundsymptomen."
        ), kind="danger")
        if rebound_t is not None
        else mo.callout(mo.md(
            "Geen rebound gedetecteerd binnen het dagvenster."
        ), kind="success")
    )
    mo.vstack([
        _callout,
        mo.accordion({
            "Hoe werkt het rebound-effect? (klik om uit te klappen)": mo.lazy(lambda: mo.md("""
    Het rebound-effect zit niet in de *hoogte* van de long tail, maar in de **snelheid van daling**.

    **Mechanisme:**
    MPH blokkeert de dopaminetransporter (DAT). Je hersenen passen zich aan aan
    verhoogde dopamine-beschikbaarheid. Wanneer de MPH-concentratie te snel daalt,
    valt de DAT-blokkade weg voordat de hersenen zijn teruggekeerd naar basislijn.
    Dit geeft tijdelijk een *relatief lagere* dopamine-beschikbaarheid dan zonder medicatie:
    prikkelbaarheid, vermoeidheid, concentratieproblemen.

    **Wanneer?**
    Rebound treedt op wanneer de curve na de dagpiek de werkzame drempel
    (hier: ~25% van piek) passeert op de dalende flank. De snelheid van die
    overgang bepaalt hoe heftig de rebound aanvoelt.

    **Wat helpt?**
    Een kleine afsluittablet (5 mg of 2.5 mg IR) voor dat moment vertraagt
    de daling en vermindert de scherpte van de overgang.

    *Bron: Volkow et al. (2003), [doi:10.1176/appi.ajp.160.11.1909](https://doi.org/10.1176/appi.ajp.160.11.1909)*
    *Innovations in CNS (2022), [doi:10.57264/cer-2022-0083](https://doi.org/10.57264/cer-2022-0083)*
            """))
        })
    ], gap="0.5rem")
    return


@app.cell(hide_code=True)
def sectie_achtergrond(mo):
    # Sectie — Achtergrond
    mo.md("""
    ## Achtergrond & literatuur
    """)
    return


@app.cell(hide_code=True)
def achtergrond(mo):
    # Achtergrond — farmacologisch model en literatuur
    mo.accordion({
        "Wat zijn IR en ER?": mo.lazy(lambda: mo.md("""
    **IR (Immediate Release) — jouw 10 mg tablet**
    Lost direct op in de maag. Piek na ~1.3 uur, werkzaam 4–5 uur.
    Curve: scherpe Bateman-piek, snel stijgend en exponentieel dalend.

    **ER (Extended Release) — jouw 27 mg Concerta**
    De buitenste laag lost direct op en geeft 22% (~6 mg) meteen vrij als IR.
    De kern pompt via osmotische druk de overige 78% (~21 mg) over ~7 uur naar buiten.
    Totale Tmax: ~6.8 uur. Curve: bifasisch (twee bulten), breed plateau.

    **AUC en dosis**
    Het oppervlak onder de curve (AUC) is proportioneel aan de totale dosis.
    25 mg IR (10+10+5) en 27 mg ER hebben vergelijkbaar oppervlak —
    het verschil is de *vorm*: smal en hoog versus breed en plat.
        """)),
        "Farmacologisch model": mo.lazy(lambda: mo.md("""
    **IR-curve:** AUC-genormaliseerde Bateman-functie.
    Parameters: ka=3.5 h⁻¹, ke=0.277 h⁻¹ → Tmax ~1.3h, t½ ~2.5h.

    **ER-curve:** Som van twee Bateman-functies:
    22% IR-coating (ka=3.5) + 78% osmotisch (ka=0.38, ke=0.198, Tmax ~3.6h).

    **Lichaamsgewicht:** ke schaalt met (W/70)^0.75 (allometrische wet).
    **LBM:** Vd schaalt met lean body mass (Janmahasatian-formule).
    **Metabolisme:** ke schaalt met exp(slider) voor CES1-variatie.
        """)),
        "Literatuur & bronnen": mo.lazy(lambda: mo.md("""
    - OROS 22%/78%: Schapperer et al. (2015), [doi:10.1002/prp2.72](https://doi.org/10.1002/prp2.72)
    - Allometrische schaling: Lyauk et al. (2016), [doi:10.1111/cts.12423](https://doi.org/10.1111/cts.12423)
    - Gewicht en klaring: Ermer et al. (2015), [doi:10.2147/DDDT.S80685](https://doi.org/10.2147/DDDT.S80685)
    - LBM formule: Janmahasatian (2005), [doi:10.1007/BF03256234](https://doi.org/10.1007/BF03256234)
    - Rebound: Volkow et al. (2003), [doi:10.1176/appi.ajp.160.11.1909](https://doi.org/10.1176/appi.ajp.160.11.1909)
    - PK-gladheid: Innovations in CNS (2022), [doi:10.57264/cer-2022-0083](https://doi.org/10.57264/cer-2022-0083)
    - Concerta SmPC: [EMA](https://www.ema.europa.eu/en/medicines/human/EPAR/concerta)
        """)),
    }, multiple=True)
    return


@app.cell(hide_code=True)
def disclaimer(mo):
    # Disclaimer
    mo.callout(mo.md(
        "**Informatief model op basis van populatie-PK parameters.** "
        "Individuele variatie door genetica (CES1), maaltijdtiming en andere medicatie "
        "kan de werkelijke curve sterk beïnvloeden. "
        "**Bespreek elk schema altijd met je arts of apotheker.**"
    ), kind="warn")
    return


if __name__ == "__main__":
    app.run()
