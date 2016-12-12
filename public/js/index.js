require('../stylesheets/index.sass')
onload=()=>{
    const trianglify = require('trianglify')
let pattern = trianglify({
    height: window.innerHeight,
    width: window.innerWidth,
    cell_size: 450,
    x_colors: ['#66CCFF', '#99CCFF', '#99FFFF'],
    y_colors: ['#339900', '#33CC99', '#33FFCC']
})
    pattern.canvas(document.getElementById('bg'))
}
let data={
    width:window.innerWidth,
    height: window.innerHeight
}