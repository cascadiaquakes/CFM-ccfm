# ccfm

[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)

Contour-based 3D fault meshing utilities. `ccfm` provides helper functions to build 3D fault meshes from contour data (top, middle and bottom LineStrings). It is part of the [CRESCENT](https://cascadiaquakes.org) Community Fault Model tools and is used by the [Contours2Mesh QGIS plugin](https://github.com/cascadiaquakes/CFM-QGIS-contours2mesh).

---

## Features

- Read/write simple contour geometries
- Geometry helpers (constants, basic operations)
- Mesh construction utilities
- Designed to be used in QGIS plugins and standalone scripts

---

## Installation

Requires Python 3.9+ (dependencies: numpy, shapely, scipy).

```bash
pip install git+https://github.com/cascadiaquakes/CFM-ccfm.git
```

---

## License

This project is licensed under the BSD 3-Clause License. See [LICENSE](LICENSE).

---

## Contacts

- **CRESCENT CFM team**: [crescentcfm@cascadiaquakes.org](mailto:crescentcfm@cascadiaquakes.org)

For bugs or questions, please [open an issue](https://github.com/cascadiaquakes/CFM-ccfm/issues).

---

## Acknowledgments

<table>
  <tr>
    <td align="center" valign="middle">
      <img src="https://raw.githubusercontent.com/cascadiaquakes/.github/main/assets/crescent-circle.png" alt="CRESCENT Logo" width="150" />
    </td>
    <td align="center" valign="middle">
      <img src="https://raw.githubusercontent.com/cascadiaquakes/.github/main/assets/nsf-logo.png" alt="NSF Logo" width="450" />
    </td>
  </tr>
</table>

This project is supported by the National Science Foundation through the Cascadia Region Earthquake Science Center (CRESCENT), Cooperative Agreement EAR-2225286.
