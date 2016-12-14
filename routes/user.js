const sql = require('./sql')

    /**
     * 
     * 
     * @param {any} req 请求头
     * @param {any} res 相应头
     */
function checklogin(req, res,fn) {
    let body = req.body
    let email=body.email
    let password=body.password
    let sqlstr='select * from user where email="'+email+'" and password="'+password+'"'
    sql.querysql(sqlstr,(err,rows)=>{
        req.session.user=rows[0]
        fn(rows[0])
    })
}
module.exports = {
    checklogin: checklogin
}