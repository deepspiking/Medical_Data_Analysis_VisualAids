# -*- coding: utf-8 -*-
"""
데이터 정정 (lossless): 연구번호 52번 mNstage 1 → 3  (정세운 선생님 확정, 2026-09-17)
========================================================================================
배경: 연구번호 52는 contra_bilateral=1(양측/대측 전이)인데 mNstage=1로 입력됨.
      Sheet1 규칙상 양측/대측 전이는 mN3 → 정정. mTstage=4이므로 mStage는
      정정 전후 모두 IV (연쇄 변화 없음).

중요: 이 xlsx에는 수식 셀(예: PFS_month 등)이 있고 openpyxl로 저장하면 수식의
      캐시 값이 소실되어 *_month 컬럼 전체가 NaN이 된다. 따라서 워크북을
      openpyxl로 재저장하지 않고, **sheet2.xml의 해당 숫자 셀만 직접 수정**한다.

동작:
  1) timestamp 백업
  2) xl/worksheets/sheet2.xml 에서 BM{row}(mNstage) 셀의 <v>1</v> → <v>3</v>
  3) 나머지 엔트리는 바이트 그대로 복사 (수식·캐시·코멘트·서식 보존)
"""
import os
import re
import shutil
import zipfile
from datetime import datetime

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "정박사님께 드릴 raw data.xlsx")
TARGET_ID = 52
SHEET_XML = "xl/worksheets/sheet1.xml"   # workbook.xml: Sheet2 -> rId1 -> sheet1.xml
COL_MN = "BM"       # mNstage (65번째 컬럼)
NEW_MN = "3"


def find_row_for_id(sheet_xml_bytes, target_id):
    """연구번호(컬럼 A)가 target_id인 XML row 번호 반환."""
    text = sheet_xml_bytes.decode("utf-8")
    rows = re.findall(r'<row[^>]*r="(\d+)"[^>]*>(.*?)</row>', text, flags=re.S)
    for rnum, body in rows:
        m = re.search(r'<c r="A' + rnum + r'"[^>]*>\s*<v>([^<]+)</v>', body)
        if m and m.group(1).strip() == str(target_id):
            return int(rnum), text
    raise SystemExit(f"연구번호 {target_id} 행을 찾지 못함")


def patch_cell_value(text, ref, new_value):
    pat = re.compile(r'(<c r="' + ref + r'"[^>]*>)(\s*<v>)([^<]*)(</v>)')
    m = pat.search(text)
    if not m:
        raise SystemExit(f"셀 {ref} 의 <v> 값을 찾지 못함")
    old = m.group(3)
    out = text[:m.start()] + m.group(1) + m.group(2) + new_value + m.group(4) + text[m.end():]
    return out, old


def main():
    bak = os.path.join(BASE, f"raw_data_backup_{datetime.now():%Y%m%d_%H%M%S}.xlsx")
    shutil.copy2(SRC, bak)
    print("백업:", os.path.basename(bak))

    with zipfile.ZipFile(SRC, "r") as zin:
        items = {n: zin.read(n) for n in zin.namelist()}

    row, text = find_row_for_id(items[SHEET_XML], TARGET_ID)
    ref = f"{COL_MN}{row}"
    text, old = patch_cell_value(text, ref, NEW_MN)
    items[SHEET_XML] = text.encode("utf-8")
    print(f"[연구번호 {TARGET_ID}] {ref} (mNstage): {old} -> {NEW_MN}")

    tmp = SRC + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for n, data in items.items():
            zout.writestr(n, data)
    os.replace(tmp, SRC)
    print("저장 완료 (lossless):", os.path.basename(SRC))


if __name__ == "__main__":
    main()
