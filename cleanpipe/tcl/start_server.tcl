proc start_server {port} {
    set server [socket -server handle_connection $port]
    puts "Server started on port $port"
    return $server
}

proc handle_connection {sock addr port} {
    puts "Connection from $addr:$port"
    fconfigure $sock -buffering line
    while {[gets $sock line] >= 0} {
        puts "Received command: $line"
        catch {eval $line} result
        puts $sock $result
        flush $sock
    }
    close $sock
}

start_server 5555
