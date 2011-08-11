#!/usr/bin/env ruby

# usage error
if ARGV.length < 1 or ARGV.length > 2 then
    $stderr.puts( "Usage: #{File.basename( $0 )} directory [cd name]" )
    exit( 1 )
end

# get variables
dir = ARGV.shift.sub( /\/$/, "" )
volname = if ARGV.length == 0 then
        File.basename( File.expand_path( dir ) )
    else
        ARGV.shift
    end
output = "#{dir}.dmg"

# dmg the folder
`hdiutil create -fs HFS+ -srcfolder "#{dir}" -volname "#{volname}" "#{output}"`
