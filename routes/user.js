const sql = require('./sql')
    /**
     * 
     * 
     * @param {any} req 请求头
     * @param {any} res 相应头
     */
function checklogin(req, res) {
    let body = req.body

    res.json(body)
}
module.exports = {
    checklogin: checklogin
}