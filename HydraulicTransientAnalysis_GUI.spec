# -*- mode: python ; coding: utf-8 -*-


block_cipher = None

from PyInstaller.utils.hooks import collect_submodules
hidden_imports_pyexcel_io = collect_submodules('pyexcel_io')
hidden_imports_numpy_core = collect_submodules('numpy')
all_hooks = hidden_imports_pyexcel_io + hidden_imports_numpy_core
a = Analysis(
    ['HydraulicTransientAnalysis_GUI.py'],
    pathex=['C:\\Users\\Sirui\\AppData\\Local\\Programs\\Python\\Python310\\lib\\site-packages'],
    binaries=[],
    datas=[],
    hiddenimports=all_hooks,
    hookspath=[r"C:\Users\Sirui\AppData\Local\Programs\Python\Python310\Lib\site-packages\numpy\_pyinstaller\hook-numpy.py"],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='HydraulicTransientAnalysis_GUI',
    debug=True,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='icon.ico',
)