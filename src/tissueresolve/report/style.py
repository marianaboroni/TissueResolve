"""
TissueResolve report design system (shared CSS).

A single, clean, publication-oriented stylesheet used by the unified report and
its components: white background, limited accent palette, strong section
headers, metric/figure/methods/warning cards, sticky sidebar navigation, and
collapsible technical details.  Import :data:`REPORT_CSS` into any page.
"""
from __future__ import annotations

__all__ = ["REPORT_CSS", "ACCENT", "ACCENT_SOFT"]

ACCENT = "#2c6e8f"        # calm teal-blue accent
ACCENT_SOFT = "#eaf2f6"

REPORT_CSS = """
:root{
  --accent:#2c6e8f; --accent-soft:#eaf2f6; --ink:#1f2933; --muted:#6b7785;
  --line:#e3e8ee; --bg:#ffffff; --card:#ffffff;
  --pass:#1e7e57; --caution:#b8860b; --warn:#c0612a; --crit:#b3261e;
  --shadow:0 1px 3px rgba(20,40,60,.07),0 1px 2px rgba(20,40,60,.04);
}
*{box-sizing:border-box}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,
  Arial,sans-serif;margin:0;color:var(--ink);background:#f6f8fa;line-height:1.55;
  font-size:15px;-webkit-font-smoothing:antialiased}
.report-container{display:flex;max-width:1180px;margin:0 auto;background:var(--bg);
  min-height:100vh;box-shadow:var(--shadow)}
/* sidebar */
.sidebar-nav{position:sticky;top:0;align-self:flex-start;height:100vh;overflow:auto;
  width:248px;min-width:248px;background:#0f2733;color:#cfe0e8;padding:20px 0;font-size:13.5px}
.sidebar-nav h1{font-size:16px;color:#fff;margin:0 18px 4px;letter-spacing:.2px}
.sidebar-nav .tagline{color:#8fb0bf;font-size:11.5px;margin:0 18px 16px}
.sidebar-nav a{display:block;color:#cfe0e8;text-decoration:none;padding:7px 18px;
  border-left:3px solid transparent}
.sidebar-nav a:hover{background:#163645;color:#fff;border-left-color:var(--accent)}
.sidebar-nav .navnum{color:#5e8597;margin-right:7px}
/* main */
main{flex:1;min-width:0;padding:28px 36px 64px}
h1.page-title{margin:0 0 2px;font-size:25px}
.page-sub{color:var(--muted);font-size:13.5px;margin:0 0 18px}
section{scroll-margin-top:16px;margin:0 0 26px}
.section-card{background:var(--card);border:1px solid var(--line);border-radius:12px;
  padding:20px 22px;box-shadow:var(--shadow)}
.section-card>h2{margin:0 0 4px;font-size:19px;border-bottom:2px solid var(--accent);
  display:inline-block;padding-bottom:3px}
h3{font-size:15px;margin:18px 0 6px;color:#243b48}
p{margin:8px 0}
a{color:var(--accent)}
.muted{color:var(--muted);font-size:12.5px}
/* metric cards */
.metric-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));
  gap:12px;margin:12px 0}
.metric-card{background:var(--accent-soft);border:1px solid var(--line);border-radius:10px;
  padding:12px 14px}
.metric-card .label{font-size:11.5px;color:var(--muted);text-transform:uppercase;
  letter-spacing:.4px}
.metric-card .value{font-size:22px;font-weight:650;color:var(--ink);margin-top:2px}
.metric-card .sub{font-size:11.5px;color:var(--muted)}
/* boxes */
.methods-card,.variable-card,.warning-card,.limitation-card,.interpretation-guide,
.how-to-read{border-radius:10px;padding:12px 15px;margin:12px 0;font-size:13.5px}
.methods-card{background:#f3f7f9;border:1px solid #d6e3ea}
.methods-card::before{content:"⚙ Methodology";display:block;font-weight:650;
  color:var(--accent);margin-bottom:4px;font-size:12px;letter-spacing:.3px}
.variable-card{background:#fbfbfd;border:1px solid var(--line)}
.variable-card::before{content:"📖 Variable definitions";display:block;font-weight:650;
  color:#54657a;margin-bottom:4px;font-size:12px}
.how-to-read{background:#fffaf0;border:1px solid #ecd9b0}
.how-to-read::before{content:"🔍 How to read this";display:block;font-weight:650;
  color:#9a7b1f;margin-bottom:4px;font-size:12px}
.interpretation-guide{background:#eef6f1;border:1px solid #cde7d8}
.interpretation-guide::before{content:"🧭 How to use this section";display:block;
  font-weight:650;color:var(--pass);margin-bottom:4px;font-size:12px}
.warning-card{background:#fff4e8;border:1px solid #e7bd95}
.warning-card::before{content:"⚠ Warnings";display:block;font-weight:650;
  color:var(--warn);margin-bottom:4px;font-size:12px}
.limitation-card{background:#f7f0f0;border:1px solid #e0c4c0}
.limitation-card::before{content:"▲ Limitations";display:block;font-weight:650;
  color:var(--crit);margin-bottom:4px;font-size:12px}
.variable-card dl{margin:4px 0;display:grid;grid-template-columns:max-content 1fr;
  gap:2px 12px}
.variable-card dt{font-weight:600;color:#2b3a47}
.variable-card dd{margin:0;color:#42505d}
/* figure cards */
.figure-card{border:1px solid var(--line);border-radius:10px;margin:14px 0;overflow:hidden}
.figure-card .fig-head{padding:10px 14px;background:#f8fafb;border-bottom:1px solid var(--line)}
.figure-card .fig-title{font-weight:650;font-size:14px}
.figure-card .fig-subtitle{color:var(--muted);font-size:12.5px}
.figure-card .fig-body{padding:6px 14px}
.figure-card img,.figure-card iframe{max-width:100%;border:0;display:block;margin:6px auto}
.figure-caption{font-size:12.5px;color:#3b4754;padding:0 14px 6px}
.figure-legend{font-size:12px;color:var(--muted);padding:0 14px 10px}
.source-data-link{display:inline-block;font-size:12px;margin:2px 10px 6px 0;
  padding:2px 9px;border:1px solid var(--line);border-radius:6px;text-decoration:none;
  background:#fafcfd}
.source-data-link:hover{background:var(--accent-soft)}
/* tables */
table{border-collapse:collapse;font-size:12.5px;margin:8px 0;width:auto}
th,td{border:1px solid var(--line);padding:4px 9px;text-align:right}
td:first-child,th:first-child{text-align:left}
th{background:#f1f5f8;font-weight:600}
.collapsible-table{margin:10px 0}
.collapsible-table>summary{cursor:pointer;font-size:13px;color:var(--accent);
  font-weight:600;padding:6px 0}
/* status badges */
.badge{display:inline-block;padding:1px 9px;border-radius:999px;font-size:11.5px;
  font-weight:650}
.status-pass{background:#e3f3ea;color:var(--pass)}
.status-caution{background:#fbf1da;color:var(--caution)}
.status-warning{background:#fbe9da;color:var(--warn)}
.status-critical{background:#f7dedb;color:var(--crit)}
.estimate-note{background:var(--accent-soft);border-left:4px solid var(--accent);
  padding:9px 13px;border-radius:6px;margin:10px 0;font-size:13px}
"""
