python -m build
if ( $LASTEXITCODE -ne 0 ) {
    Write-Error -Message "Build failed" -ErrorAction Stop
}

twine upload -r local dist/*
cd ../rotsim2d_apps
pip install --extra-index-url=http://127.0.0.1:4040/simple --no-deps --pre --upgrade --force-reinstall --no-cache-dir rotsim2d molspecutils
cd ../rotsim2d
