# -*- coding: utf-8 -*-
"""
nomogram 산출 그림을 논문용 패널(A/B/C)로 합치기
출력: results/exp1/nomogram/panels/
"""
import os
from PIL import Image, ImageDraw

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "results", "exp1", "nomogram")
OUT = os.path.join(SRC, "panels")
os.makedirs(OUT, exist_ok=True)

PAD = 24
LABEL_H = 46


def hstack(paths, labels, out, bg="white"):
    imgs = [Image.open(p).convert("RGB") for p in paths]
    h = max(i.height for i in imgs)
    imgs = [i.resize((int(i.width * h / i.height), h), Image.LANCZOS) for i in imgs]
    W = sum(i.width for i in imgs) + PAD * (len(imgs) + 1)
    H = h + PAD * 2 + LABEL_H
    canvas = Image.new("RGB", (W, H), bg)
    draw = ImageDraw.Draw(canvas)
    x = PAD
    for img, lab in zip(imgs, labels):
        canvas.paste(img, (x, PAD + LABEL_H))
        draw.text((x + 6, PAD + 6), lab, fill="black")
        x += img.width + PAD
    canvas.save(out, dpi=(300, 300))
    print("[save]", os.path.basename(out), canvas.size)


for y in ["OS", "PFS"]:
    hstack([os.path.join(SRC, f"calibration_{y}_3yr.png"),
            os.path.join(SRC, f"calibration_{y}_5yr.png")],
           ["A. 3-year calibration", "B. 5-year calibration"],
           os.path.join(OUT, f"panel_{y}_calibration.png"))
    hstack([os.path.join(SRC, f"dca_{y}_3yr.png"),
            os.path.join(SRC, f"roc_{y}.png"),
            os.path.join(SRC, f"risk_km_{y}.png")],
           [f"A. Decision curve (3-year)", "B. Time-dependent ROC",
            "C. Risk-group KM"],
           os.path.join(OUT, f"panel_{y}_decision.png"))
print("[출력]", OUT)
