# rev27 — REL_BENT (offset bend) 세그먼트 추가 요약

## 변경점 (rev26 → rev27)

| 위치 | 추가/수정 |
|---|---|
| 헤더 docstring | 버전 rev27, 세그먼트 목록에 `REL_BENT` 추가, rev27 기능 소개 |
| `Segment` dataclass 주석 | kind 목록에 `'rel_bent'` 추가 |
| `segment_geometry()` | `kind == "rel_bent"` geometry 블록 신설 (두 원호 합성) |
| `SEG_PARAMS["rel_bent"]` | `radius`, `dy` 두 파라미터 등록 (Add Segment 폼 노출) |
| `SEG_DISPLAY["rel_bent"]` | 드롭다운 표시명 추가 |
| `_segment_signed_curvature()` | straight-like 그룹에 `rel_bent` 추가 (양 끝 κ=0 취급) |
| `_show_about()` | About 다이얼로그에 rev27 / REL_BENT 표기 |

## REL_BENT 란?

굽힘 반경 **R** 과 측면 오프셋 **dy** 만 주면, **동일한 두 개의 원호**로
S-bend(offset bend)를 만들어 **끝단 진행각이 항상 시작각(0°)으로 복귀**하는
세그먼트입니다.

- 첫 번째 호: 곡률 `+s/R`, `+s·θ` 회전
- 두 번째 호: 곡률 `−s/R`, `−s·θ` 회전 → 시작 heading 으로 복귀
- `s = sign(dy − y0)`

즉 시작 포트와 동일한 방향을 보면서 **절대 y 좌표 `dy`** 에 도달하는
끝 포트를 자동 생성합니다. SBEND 와 달리 곡률 반경이 고정값으로
직접 지정되므로, 최소 굽힘 반경 제약(곡률 손실)을 보장하기 좋습니다.

## 절대 y 타깃 + `yoffset(n)` (rev27 업데이트)

`dy` 는 **상대 변위가 아니라 도달할 절대 y 좌표**로 해석됩니다.
실제 적용되는 측면 변위는

```
h_signed = dy − y0        (y0 = 유효 시작 y)
```

따라서 `dy = yoffset(5)` 라고 쓰면 **segment 5 의 끝점 절대 y** 에
정확히 정렬되는 bend 가 만들어집니다. `yoffset(n)` 은 모든 숫자
입력 필드/파라미터 식에서 쓸 수 있는 헬퍼로, segment n 의 끝점
포트 y 를 돌려줍니다(참조 대상 segment 가 시퀀스상 먼저 있어야 함).
절대-y 타깃은 시작 heading 이 대체로 수평(0°)인 일반적 사용을 가정합니다.

## 물리 / 기하 공식

```
h = |dy − y0|              (적용할 측면 변위)
h = 2 R (1 − cos θ)        →     θ = arccos(1 − h / (2R))
전진거리  dx = 2 R sin θ      (R, h 로 자동 결정)
```

- 각 호의 회전각 **θ** 는 R 과 h 로 완전히 결정됩니다.
- 도달 가능 범위: `h ≤ 4R` (θ ≤ 180°). 초과 시 `h = 4R` 로
  클램프되어 θ = 180° 가 됩니다.
- `h = 0`(dy == y0) 또는 `R ≤ 0` 이면 straight stub(퇴화) 으로 처리됩니다.

## 파라미터

| key | 라벨 | 기본값 | 의미 |
|---|---|---|---|
| `radius` | Bend radius [µm] | 500.0 | 두 원호의 굽힘 반경 R |
| `dy` | Target y [µm] (abs; e.g. yoffset(5)) | 100.0 | 도달할 **절대 y 좌표** |

(`width` / `xoffset` / `yoffset` / `angleoffset` 등 공통 오프셋 파라미터도
다른 세그먼트와 동일하게 적용됨. 단, rel_bent 의 `dy` 는 절대 y 타깃이므로
공통 `yoffset` 파라미터로 시작점을 옮기면 그만큼 h 가 재계산됨)

## Auto-Offset 와의 관계

`rel_bent` 는 시작·끝 heading 이 동일(0° 복귀)하므로 `sbend` 와 동일하게
**straight-like** 로 취급됩니다(`_segment_signed_curvature` 가 양 끝단 κ=0
반환). 따라서 rev25 의 auto bend-offset 접합 보정에서 추가 offset 을
유발하지 않습니다.

## 검증

- 끝점 y 가 시작 y0 와 무관하게 목표 **절대 y `dy`** 와 정확히 일치
  (예: y0=50, dy=200 → end_y=200.000 / y0=30, dy=500 → end_y=500.000).
- `yoffset(n)`: segment n 끝점 y 반환, 식 `dy="yoffset(5)"` 로 정렬 동작,
  미존재 segment 는 명확한 에러.
- 클램프 케이스(h > 4R) → h = 4R, θ = 180° 로 정상 제한.
- 끝 heading 은 항상 시작각(`end_a = a0`).
- `python3 -m py_compile` 구문 검사 통과.

## 하위 호환

rev26 이전 프로젝트 파일(JSON)에는 `rel_bent` 세그먼트가 없으므로
영향이 없으며, 기존 세그먼트 종류·동작은 완전히 보존됩니다.
